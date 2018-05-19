#ifndef MPC_SOLVER_HPP
#define MPC_SOLVER_HPP

#include "mpc_solver.h"
#include "Timer.h"
//#include "Utilities/Math/PseudoInverse.h"

template<int M, int N>
MPC_Solver<M,N>::MPC_Solver(double integration_tol, double solver_tol, double solTime):
    //nr_solver( this, &MPC_Solver::costJacobian, &MPC_Solver::costHessian ),
    p_trajectory_definition(0)
{
    validation_mode = false;

    solveTime = solTime;

    int_tol = integration_tol;
    nr_tol = solver_tol;
}

template<int M, int N>
bool MPC_Solver<M,N>::solve(const Eigen::Matrix<double,M,1>& Xo_, std::vector< Trajectory<M,N> >& trajectory_definition, const std::vector< Dynamic_System_Constraint<M,N>* >& constraint_list, const std::vector< Dynamic_System_Constraint<M,N>* >& soft_constraint_list_)
{
    //    Timing::stopWatch timer;
    typename std::vector< Trajectory<M,N> >::iterator it; // trajectory iterator

    this->constraint_list = constraint_list;
    this->soft_constraint_list = soft_constraint_list_;



    X_IC = Xo_;
    // save trajectory definition pointer for use in other internal functions
    p_trajectory_definition = &trajectory_definition;

    numConstraints  = trajectory_definition.size() * this->constraint_list.size();
    numControls = trajectory_definition.size() *N; // no control for the IC
    numStates = trajectory_definition.size() *M; // no control for the IC

    // Initialize System Control inputs to be target controls
    bool resolve = false;
    U.resize(numControls+numConstraints,1);
    U_Target.setZero(numControls,1);
    Ju.setZero(numStates, numControls);
    Jx.setIdentity(numStates, numStates);


    nontrivial_Yu = false;
    nontrivial_Yxx = false;
    nontrivial_S = false;
    nontrivial_Yuu = false;
    nontrivial_Yux = false;
    nontrivial_H = false;
    nontrivial_Cux = false;

    numOutputs = 0;
    it = trajectory_definition.begin();
    for(int rowcountTarget=0; it!= p_trajectory_definition->end(); it++,rowcountTarget+=N)
    {
        if( it->includeTargetControlInCost )//for kalman filtering
        {
            U_Target.block(rowcountTarget,0,N,1) = it->target_control;
        }

        if( it->soln_control.rows() == N)
        {
            U.block(rowcountTarget,0,N,1) = it->soln_control;
            resolve = true;
        }
        else
            U.block(rowcountTarget,0,N,1) = 0*Eigen::Matrix<double,N,1>::Ones();

        numOutputs += it->p_dynamic_system->getNumOutputs();
        nontrivial_Yu |= it->p_dynamic_system->has_out_u_sense();
        nontrivial_Yxx |= it->p_dynamic_system->has_out_xx_sense();
        nontrivial_S |= it->p_dynamic_system->has_out_xx_sense() || it->p_dynamic_system->has_out_ux_sense() || it->p_dynamic_system->has_ode_xu_sense() || it->p_dynamic_system->has_ode_uu_sense();
        nontrivial_Yuu |= it->p_dynamic_system->has_out_uu_sense();
        nontrivial_Yux |= it->p_dynamic_system->has_out_ux_sense();
        nontrivial_H |= it->p_dynamic_system->has_out_xx_sense() || it->p_dynamic_system->has_ode_xx_sense();
    }


    if( this->constraint_list.size() )
    {
        //typename std::vector< Dynamic_System_Constraint<M,N>* >::iterator it_c; // trajectory iterator

        for(unsigned int cons_ind = 0; cons_ind < this->constraint_list.size(); cons_ind++)
        {
            Dynamic_System_Constraint<M,N>* dsc = this->constraint_list[cons_ind];
            nontrivial_H |= dsc->has_xx_sens();
            nontrivial_Cux |= dsc->has_ux_sens();
        }

        if( !resolve )
        {
            Eigen::VectorXd ones; ones.setOnes(numConstraints,1);
            U.block(numControls,0,numConstraints,1) = 0*ones;
        }
    }

    if( soft_constraint_list.size() )
    {
        for(unsigned int cons_ind = 0; cons_ind < this->soft_constraint_list.size(); cons_ind++)
        {
            Dynamic_System_Constraint<M,N>* dsc = this->soft_constraint_list[cons_ind];
            nontrivial_H |= dsc->has_xx_sens();
            nontrivial_Cux |= dsc->has_ux_sens();
        }

    }

    //Set Q and R matrices to zeros
    Q.setZero(numOutputs,numOutputs);
    R.setZero(numControls,numControls);


    Jg.setZero(numOutputs, numStates);
    E.setZero(1,numOutputs);

    Z_constraint.setZero(1,numStates);
    Z.setZero(1,numStates);            //Corrected from being numOutput to numState
    Y.setZero(1,numOutputs);

    //Packing Q
    it = p_trajectory_definition->begin();
    for(int diagElement=0; it!= p_trajectory_definition->end(); it++)
    {
        int stepSumOutput = it->p_dynamic_system->getNumOutputs();
        Q.block(diagElement,diagElement,stepSumOutput,stepSumOutput) = it->Q;
        diagElement+=it->p_dynamic_system->getNumOutputs();
    }
    //cout << "Q: " << endl << Q << endl << endl;

    //Packing R
    it = p_trajectory_definition->begin();
    for(int diagElement=0; it!= p_trajectory_definition->end(); it++,diagElement+=N)
    {
        R.block(diagElement,diagElement,N,N) = it->R;
    }

    /*
     //Packing C
     C.setZero(numOutputs,numStates);

     int rowNum = 0;
     int colNum = 0;
     for( it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
     {
         int nOutput = it->p_dynamic_system->getNumOutputs();
         C.block(rowNum,colNum,nOutput,M) = it->p_dynamic_system->getOutSelectMatrix();
         rowNum += nOutput;
         colNum += M;
     }

     // pack CtQC;
     CtQC = C.transpose()*Q*C;
     */


    H.setZero(numStates,numStates);
    S.setZero(numControls,numControls);
    Joo.setZero(numControls,numControls);
    Jox.setZero(numControls,numStates);

    // solve system
    //    cout << "Initializtion Time: " << timer.stop()<<endl;
    //nr_solver.solve(U, tolerance, 1000,false,NEWTON_RAPHSON,1e-4);

    //    if(  !resolve  )
    //    {
    //        this->constraint_list.clear();
    //        //nr_solver.solve(U, nr_tol, 500,true,LEVENBERG_MARQUARDT_ADAPTIVE,1);
    //        newton_raphson();
    //        this->constraint_list = constraint_list;
    //    }


    bool converged = newton_raphson();
    //    cout << "Solve Time: " << timer.stop()<<endl;
    //nr_solver.solve(U,1e-9,1000,true,LEVENBERG_MARQUARDT_ADAPTIVE,1e-4);

    // extract solution from solver
    //  U = nr_solver.solution();

    // unpack solution
    it = trajectory_definition.begin();
    for(int rowElement=0; it!= p_trajectory_definition->end(); it++,rowElement+=N)
    {
        it->soln_control = U.block(rowElement,0,N,1) + it->target_control;//U_Target.block(rowElement,0,N,1);
        it->soln_state = it->p_dynamic_system->x();
        it->soln_output = it->p_dynamic_system->getOutSelectMatrix()*it->p_dynamic_system->y();
    }

   // cout << "First Step Target Control: " << p_trajectory_definition->front().target_control.transpose() << endl;
    //cout << "First Step Control: " << p_trajectory_definition->front().soln_control.transpose() << endl;
    // save trajectory definition pointer to zero because it is no longer valid to this class
    p_trajectory_definition = 0;
    this->constraint_list.clear();

    //    cout << "Cleanup Time: " << timer.stop()<<endl;

    return converged;

}

template<int M, int N>
double MPC_Solver<M,N>::cost(  const Eigen::VectorXd& current_controls )
{
    double ret_cost = 0;

    typename std::vector< Trajectory<M,N> >::iterator it; // trajectory iterator

    Eigen::Matrix<double,M,1> Xo = X_IC;
    int j=0;
    int controlInd = 0;
    int lambdaIndex = numControls;
    for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++ )
    {
        //Computing my intial y here.
        it->p_dynamic_system->solve(Xo, current_controls.block(controlInd,0,N,1) + it->target_control, it->t_span, true, 1e-8, 1000 );
        //int segNumOutput = it->p_dynamic_system->getNumOutputs();

        Eigen::VectorXd curControl = current_controls.block(controlInd,0,N,1) + U_Target.block(controlInd,0,N,1);
        //if( it->includeTargetControlInCost )
        //    curControl += it->target_control;

        Eigen::VectorXd outputError = it->p_dynamic_system->getOutSelectMatrix()*it->p_dynamic_system->y() - it->target_output;

        Eigen::MatrixXd tmpCost = (outputError.transpose()* it->Q*outputError + curControl.transpose()*it->R*curControl);
        ret_cost += tmpCost(0,0);

        Xo = it->p_dynamic_system->x();// update based on last solution....


        if( this->constraint_list.size() )
        {

            for( unsigned int cons_ind = 0; cons_ind < constraint_list.size(); cons_ind++)
            {
                Dynamic_System_Constraint<M,N>* dsc = constraint_list[cons_ind];
                ret_cost += current_controls(lambdaIndex)*dsc->constraint_cost(it->t_span.tf,Xo,curControl, HARD);
                lambdaIndex++;
            }
        }

        if( this->soft_constraint_list.size() )
        {

            for( unsigned int cons_ind = 0; cons_ind < soft_constraint_list.size(); cons_ind++)
            {
                Dynamic_System_Constraint<M,N>* dsc = soft_constraint_list[cons_ind];
                ret_cost += dsc->constraint_cost(it->t_span.tf,Xo,curControl, SOFT);
            }
        }

        controlInd += N;
    }

    return ret_cost;
}

template<int M, int N>
void MPC_Solver<M,N>::costJacobian( Eigen::VectorXd& current_controls, Eigen::VectorXd& return_costJacobian)
{
    //   Timing::stopWatch timer;
    assert(p_trajectory_definition && "MPC_Solver::costJacobian Error: Trajectory pointer not set");


    // resize costJacobian as requried
    return_costJacobian.setZero(numControls+numConstraints,1);

    typename std::vector< Trajectory<M,N> >::iterator it; // trajectory iterator

    int rowNum = 0;
    int colNum = 0;

    it = p_trajectory_definition->begin();
    U = current_controls;

    Eigen::Matrix<double,M,1> Xo =X_IC;
    int j=0;
    for(int i=0; it!= p_trajectory_definition->end(); it++,i+=N )
    {
        //Computing my intial y here.
        it->p_dynamic_system->solve(Xo, U.block(i,0,N,1) + it->target_control, it->t_span, true, int_tol, 10000 );
        int segNumOutput = it->p_dynamic_system->getNumOutputs();

        // calc error
        E.block(0,j,1,segNumOutput) = ((it->p_dynamic_system->getOutSelectMatrix() * it->p_dynamic_system->y()) - it->target_output).transpose();

        Xo = it->p_dynamic_system->x();// update based on last solution....
        j += segNumOutput;
    }
    //    cout << "Integration Time (more or less): " << timer.stop()<<endl;

    //   cout << "E: " << E << endl;

    //Packing Ju

    int contCount = 0;
    int numRow=0;
    for( it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
    {
        Ju.block(numRow,contCount,M,N) = it->p_dynamic_system->Ju_sens();
        contCount += N;
        numRow += M;
    }
    /*contCount = N;
    for(int numRow = M; numRow < numStates; numRow += M )
    {
        Ju.block(numRow,contCount,(numStates-numRow),N) = Ju.block(numRow,contCount-N,(numStates-numRow),N);
        contCount += N;
    }*/

    //Packing Jx
    ///We need +1 here for trajectory begin so that we start with Jx2 right?
    it = p_trajectory_definition->begin()+1; // start with Jx2 not Jx1
    for(int rowCount=M; it!= p_trajectory_definition->end(); it++, rowCount+=M )
    {
        for (int colmCount = (rowCount-M); colmCount >=0; colmCount -= M)   // only access lower triangle
        {
            if( rowCount-M == colmCount )
            {
                // identity multiply
                Jx.block(rowCount,colmCount, M, M) = it->p_dynamic_system->Jx_sens();
            }
            else
            {
                Jx.block(rowCount,colmCount, M, M) = it->p_dynamic_system->Jx_sens() * Jx.block(rowCount-M,colmCount, M, M);
            }
        }
    }

    //cout << "Jx" << endl << Jx << endl;

    //Packing Jg
    rowNum = 0;
    colNum = 0;
    for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
    {
        int numOutput = it->p_dynamic_system->getNumOutputs();

        Jg.block(rowNum, colNum, numOutput, M) = it->p_dynamic_system->getOutSelectMatrix()*it->p_dynamic_system->Yx_sens();

        colNum += M;
        rowNum += numOutput;
    }

    //cout << Jg << endl;



    Y = E*Q;
    // cout << "Y: " << Y << endl;
    Z = Y*Jg*Jx;
    // cout << "Z: " << Z << endl;

    JxJu = Jx*Ju;

    //cout << "JgJxJu:\n"<< Jg*JxJu << endl <<endl;

    return_costJacobian.block(0,0,numControls,1) = 2.0*((Z*Ju).transpose() + R*(U.block(0,0,numControls,1)+U_Target));

    // cout << "cost: " << return_costJacobian.transpose() << endl;
    // cout << endl << "Ju:"<<Ju << endl << endl;
    //Packing Jo
    if( nontrivial_Yu )
    {
        rowNum=0;
        for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
        {
            int numOutput = it->p_dynamic_system->getNumOutputs();
            Jo.block(rowNum, 0, numOutput, N) = it->p_dynamic_system->getOutSelectMatrix()*it->p_dynamic_system->Yu_sens();
            rowNum += numOutput;
        }
        return_costJacobian.block(0,0,numControls,1) += 2.0*(Y*Jo).transpose();
    }

    // add in cost constraints
    if( this->constraint_list.size() )
    {
        Z_constraint.setZero(1,numStates);
        x_der_constraints.setZero(numConstraints,numStates);
        u_der_constraints.setZero(numConstraints,numControls);
        int costIndex = numControls;
        int trajIndex = 0;

        for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
        {
            double time = it->t_span.tf;
            Eigen::Matrix<double, M, 1> state_current( it->p_dynamic_system->x() );
            Eigen::Matrix<double, N, 1> control_current( U.block(trajIndex*N,0,N,1) + it->target_control );
            for( unsigned int cons_ind = 0; cons_ind < constraint_list.size(); cons_ind ++ )
            {
                Dynamic_System_Constraint<M,N>* dsc = constraint_list[cons_ind];
                return_costJacobian(costIndex) = dsc->constraint_cost(time, state_current, control_current,  HARD);//constraint cost driven to zero by tf.

                if( dsc->has_u_sens() )
                {
                    u_der_constraints.block(costIndex-numControls,trajIndex*N,1,N) = dsc->u_sens(time, state_current, control_current,  HARD).transpose();
                    // cout << u_der_constraints << endl << endl;
                    return_costJacobian.block(trajIndex*N,0,N,1) += U(costIndex)*u_der_constraints.block(costIndex-numControls,trajIndex*N,1,N).transpose();

                }
                if( dsc->has_x_sens() )
                {
                    x_der_constraints.block(costIndex-numControls,trajIndex*M,1,M) = dsc->x_sens(time, state_current, control_current,  HARD).transpose();
                    Z_constraint.block(0,trajIndex*M,1,M) += U(costIndex)*x_der_constraints.block(costIndex-numControls,trajIndex*M,1,M);
                }


                costIndex ++;
            }
            trajIndex ++;
        }

        Z_constraint = Z_constraint*Jx;
        return_costJacobian.block(0,0,numControls,1) += (Z_constraint*Ju).transpose();
    }

    if( soft_constraint_list.size() )
    {
        Z_soft_constraint.setZero(1,numStates);
        int trajIndex = 0;

        for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
        {
            double time = it->t_span.tf;
            Eigen::Matrix<double, M, 1> state_current( it->p_dynamic_system->x() );
            Eigen::Matrix<double, N, 1> control_current( U.block(trajIndex*N,0,N,1) + it->target_control );
            for( unsigned int cons_ind = 0; cons_ind < soft_constraint_list.size(); cons_ind ++ )
            {
                Dynamic_System_Constraint<M,N>* dsc = soft_constraint_list[cons_ind];
                if( dsc->has_u_sens() )
                {
                    // cout << u_der_constraints << endl << endl;
                    return_costJacobian.block(trajIndex*N,0,N,1) += dsc->u_sens(time, state_current, control_current,  SOFT);

                }
                if( dsc->has_x_sens() )
                {
                    Z_soft_constraint.block(0,trajIndex*M,1,M) += dsc->x_sens(time, state_current, control_current,  SOFT).transpose();
                }

            }
            trajIndex ++;
        }

        Z_soft_constraint = Z_soft_constraint*Jx;
        return_costJacobian.block(0,0,numControls,1) += (Z_soft_constraint*Ju).transpose();
    }

    //    if( !validation_mode )
    //    {
    //        cout <<"cost: " << return_costJacobian.transpose() << endl;
    //        cout << "current controls: " << current_controls.transpose() << endl;
    //    }

    //    cout << "Jacobian Calculation Time: " << timer.stop()<<endl;
    return;
}
template<int M, int N>
void MPC_Solver<M,N>::costHessian( Eigen::VectorXd& current_controls, Eigen::MatrixXd& return_costHessian )
{
    //    Timing::stopWatch timer;

    // verify costJacobian was called first
    assert( (U - current_controls).norm() == 0 && "Current Control Vector Not as Expected" );


    // zero out and resize costHessian if necessary
    return_costHessian.setZero(numControls + numConstraints , numControls + numConstraints );

    typename std::vector< Trajectory<M,N> >::iterator it; // trajectory iterator
    typename std::vector< Dynamic_System_Constraint<M,N>* >::iterator it_c; // trajectory iterator

    if( constraint_list.size() )
        Z += Z_constraint/2.0;

    if( soft_constraint_list.size() )
        Z += Z_soft_constraint/2.0;

    if( nontrivial_Cux )
        Cux.setZero(numControls,numStates);


    //Pack H Matrix
    if( nontrivial_H )
    {
        int diagNum = 0;                   //We need from Z(i+1)
        Eigen::RowVectorXd zBlock;  //defined separately to help the compiler out
        for(it = p_trajectory_definition->begin()+1; it!= p_trajectory_definition->end(); it++)    //Since we need Jxx(i+1)
        {
            if( it->p_dynamic_system->has_ode_xx_sense() )
            {
                zBlock = Z.block(0,diagNum+M,1,M);
                H.block(diagNum,diagNum,M,M) = zBlock * it->p_dynamic_system->Jxx_sens();

                if( it->p_dynamic_system->has_out_xx_sense() )
                {
                    zBlock = Y.block(0,diagNum,1,M)*it->p_dynamic_system->getOutSelectMatrix();
                    H.block(diagNum,diagNum,M,M) += zBlock * it->p_dynamic_system->Yxx_sens();
                }

            }
            else if( it->p_dynamic_system->has_out_xx_sense() )
            {
                zBlock = Y.block(0,diagNum,1,M)*it->p_dynamic_system->getOutSelectMatrix();
                H.block(diagNum,diagNum,M,M) = zBlock * it->p_dynamic_system->Yxx_sens();
            }


            diagNum += M;
        }
    }


    //Pack S Matrix
    if( nontrivial_S )
    {
        S.setZero(numControls,numControls);
        int z_rowNum = 0;
        int s_rowNum = 0;
        int colNum = 0;

        Eigen::RowVectorXd zBlock;  //defined separately to help the compiler out
        for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
        {
            zBlock = Z.block(0,z_rowNum,1,M);
            for(colNum = 0; colNum <= s_rowNum; colNum++)
            {
                if (colNum==s_rowNum && it->p_dynamic_system->has_ode_uu_sense())        //Juu sensitivity terms on the diagonal
                {
                    S.block(s_rowNum,colNum,N,N) += 2.0*zBlock * it->p_dynamic_system->Juu_sens();
                    //cout << endl << S.block(s_rowNum,colNum,N,N) << endl << endl;
                }

                if( it->p_dynamic_system->has_ode_xu_sense() && z_rowNum > 0) // >M?
                {
                //    cout << it->p_dynamic_system->Jxu_sens() << endl <<endl;
                    Eigen::MatrixXd JxuJu = (zBlock * it->p_dynamic_system->Jxu_sens())* JxJu.block(z_rowNum-M,colNum,M,N);
                    S.block(s_rowNum, colNum,N,N) += JxuJu;
                    S.block(colNum, s_rowNum,N,N) += JxuJu.transpose();
                //    cout << S.block(s_rowNum, colNum,N,N) << endl << endl;
                }
            }
            s_rowNum+=N;
            z_rowNum+=M;
        }
    }
   // cout << "Z: " << Z << endl;
   // cout << "S:\n" << S << endl <<endl;


    if( nontrivial_Yu && (nontrivial_Yuu || nontrivial_Yux) )
    {
        //Pack Joo and Jox
        int rowNum  = 0;
        int diagNum = 0;
        int colNum  = 0;
        for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
        {
            int nOut = it->p_dynamic_system->getNumOutputs();
            Eigen::RowVectorXd yBlock(Y.block(0,rowNum,1,nOut) * it->p_dynamic_system->getOutSelectMatrix());

            if( it->p_dynamic_system->has_out_uu_sense() )
            {
                Joo.block(diagNum, diagNum, N,N) = yBlock * it->p_dynamic_system->Yuu_sens();     //How do we include the selection matrix multiplication?
            }

            if( it->p_dynamic_system->has_out_ux_sense() )
            {
                Jox.block(rowNum,colNum,M,nOut)  = yBlock * it->p_dynamic_system->Yux_sens();
            }

            diagNum += N;
            rowNum += nOut;
        }

    }



    //total Hessian so far..
    //JxJu = Jx*Ju; Moved to cost size for constraints

    // add in cost constraints
    if( this->constraint_list.size() )
    {
        int costIndex = numControls;
        int trajIndex = 0;
        int s_rowNum = 0;

        for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
        {
            double time = it->t_span.tf;
            Eigen::Matrix<double, M, 1> state_current( it->p_dynamic_system->x() );
            Eigen::Matrix<double, N, 1> control_current( U.block(trajIndex*N,0,N,1) + it->target_control );
            for( unsigned int cons_ind = 0; cons_ind < constraint_list.size(); cons_ind ++ )
            {
                Dynamic_System_Constraint<M,N>* dsc = constraint_list[cons_ind];
                if( dsc->has_xx_sens() ) //0.5* is because H will be multiplied by 2 later
                    H.block(trajIndex*M,trajIndex*M,M,M) += 0.5*U(costIndex)*dsc->xx_sens(time,state_current,control_current,  HARD);

                if( dsc->has_uu_sens())
                    return_costHessian.block(trajIndex*N,trajIndex*N,N,N) += U(costIndex)*dsc->uu_sens(time,state_current,control_current,  HARD);


                if( dsc->has_ux_sens())
                    Cux.block(trajIndex*N,trajIndex*M,N,M) += U(costIndex)*dsc->ux_sens(time,state_current,control_current,  HARD);

                costIndex ++;
            }
            trajIndex ++;
            s_rowNum+=N;
        }

        // implement symmetry
        return_costHessian.block(numControls,0,numConstraints,numControls) = x_der_constraints*JxJu + u_der_constraints;
        // cout << endl <<  return_costHessian.block(numControls,0,numConstraints,numControls) << endl;
        return_costHessian.block(0,numControls,numControls,numConstraints) = return_costHessian.block(numControls,0,numConstraints,numControls).transpose();
    }

    if( soft_constraint_list.size() )
    {
        int trajIndex = 0;
        int s_rowNum = 0;

        for(it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it++)
        {
            double time = it->t_span.tf;
            Eigen::Matrix<double, M, 1> state_current( it->p_dynamic_system->x() );
            Eigen::Matrix<double, N, 1> control_current( U.block(trajIndex*N,0,N,1) + it->target_control );
            for( unsigned int cons_ind = 0; cons_ind < soft_constraint_list.size(); cons_ind ++ )
            {
                Dynamic_System_Constraint<M,N>* dsc = soft_constraint_list[cons_ind];
                if( dsc->has_xx_sens() ) //0.5* is because H will be multiplied by 2 later
                    H.block(trajIndex*M,trajIndex*M,M,M) += 0.5*dsc->xx_sens(time,state_current,control_current,  SOFT);

                if( dsc->has_uu_sens())
                    return_costHessian.block(trajIndex*N,trajIndex*N,N,N) += dsc->uu_sens(time,state_current,control_current,  SOFT);


                if( dsc->has_ux_sens())
                    Cux.block(trajIndex*N,trajIndex*M,N,M) += dsc->ux_sens(time,state_current,control_current,  SOFT);


            }
            trajIndex ++;
            s_rowNum+=N;
        }

    }

    if(nontrivial_Cux)
    {
        Eigen::MatrixXd CuxJxJu(Cux*JxJu);
        return_costHessian.block(0,0,numControls,numControls) += CuxJxJu + CuxJxJu.transpose();
    }

    if( nontrivial_H )
    {
        return_costHessian.block(0,0,numControls,numControls) +=  2.0*(R + (JxJu).transpose()*( Jg.transpose()*Q*Jg + H)*JxJu);
   //     cout << "--\n"<<return_costHessian.block(0,0,numControls,numControls) << endl<<endl;
    }
    else
    {
        JgJxJu = Jg*JxJu;
        return_costHessian.block(0,0,numControls,numControls) +=  2.0*(R + JgJxJu.transpose()*Q*JgJxJu);
    }

    if ( nontrivial_S )
    {
        return_costHessian.block(0,0,numControls,numControls) += S;
   //     cout << "--\n"<<return_costHessian.block(0,0,numControls,numControls) << endl<<endl;
    }

    if ( nontrivial_Yu )
        return_costHessian.block(0,0,numControls,numControls) += 2.0*Jo.transpose()*Q*Jo;

    if( nontrivial_Yuu )
        return_costHessian.block(0,0,numControls,numControls) += 2.0*Joo;

    if( nontrivial_Yux )
    {
        JoxJxJu = Jox*JxJu;
        return_costHessian.block(0,0,numControls,numControls) += JoxJxJu + JoxJxJu.transpose();
    }

    //    if( !validation_mode )
    //        validate_jacobians( X_IC, U, *p_trajectory_definition, constraint_list );

    //    cout << "Packing Hessian Time: " << timer.stop() <<endl;
    return;
}

template<int M, int N>
bool MPC_Solver<M,N>::newton_raphson()
{

    // nr_solver.solve(U, nr_tol, 1000,true,LEVENBERG_MARQUARDT_ADAPTIVE,1e-4);
    // nr_solver.solve(U, nr_tol, 1000,true,LEVENBERG_MARQUARDT_ADAPTIVE,1e-6);
    solveTimer.start();

    Eigen::MatrixXd Hessian;
    Eigen::MatrixXd HessianReduced;
    Eigen::MatrixXd Scale;


    Scale.setZero(numControls + numConstraints , numControls + numConstraints );

    Eigen::VectorXd Error;
    Eigen::VectorXd ErrorReduced;
    Eigen::VectorXd dU;

    Eigen::VectorXd U_last = U;

    std::vector<unsigned int> nonZeroCols;
    std::vector<double> colNorms;

    Eigen::LDLT<Eigen::MatrixXd> Solver;
    //Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Solver;
    //Eigen::FullPivLU<Eigen::MatrixXd> Solver;
    //Eigen::LLT<Eigen::MatrixXd> Solver;

    //Solver.setThreshold(1e-9);

    U.tail(numConstraints) *= 0;
    double Lambda = 1e-6;
    costJacobian(U,Error);

    int iterations = 0;

    do
    {
        nonZeroCols.clear();
        costHessian(U,Hessian);

//        double solveTime = solveTimer.lap();

//        Eigen::VectorXd uTmp(U.head(numControls+p_trajectory_definition->size()*constraint_list.size()));
//        if( !validate_jacobians( X_IC, uTmp, *p_trajectory_definition, constraint_list, soft_constraint_list,true ))
//        {
//            for( unsigned int i=0; i<p_trajectory_definition->size(); i++ )
//                (*p_trajectory_definition)[i].p_dynamic_system->verify_jacobians();
//        }

//        if( constraint_list.size() && !constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N),HARD) )
//            int i=0;
//        if( soft_constraint_list.size() && !soft_constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N),SOFT) )
//            int i=0;

        //  cout << "Hessian Orig:\n"<<Hessian<<endl<<endl;



        double maxNorm = 0;
        colNorms.clear();
        for( unsigned int i=0; i<Hessian.cols(); i++)
        {
            colNorms.push_back(Hessian.col(i).norm());
            maxNorm = std::max(colNorms.back(),maxNorm);
        }
        for( unsigned int i=0; i<Hessian.cols(); i++)
        {
            if(colNorms[i] > 0)
            {
                nonZeroCols.push_back(i);
            }
        }
        if( nonZeroCols.size() == 0 )
            return 0;

        Eigen::MatrixXd reducer;

        reducer.setZero(Hessian.rows(),nonZeroCols.size());
        int j=0;
        for( unsigned int i=0; i<Hessian.cols() && j<nonZeroCols.size(); i++)
        {
            if( nonZeroCols[j] == i )
            {
                reducer(i,j) = 1;
                j++;
            }
        }



        //  dU = Solver.compute(scale*reducer.transpose()*Hessian*reducer).solve(scale*reducer.transpose()*Error);
        //  double threshold = Solver.threshold();

        HessianReduced = reducer.transpose()*Hessian*reducer;


       // Scale.setZero(HessianReduced.rows(),HessianReduced.cols());

        for( int i = 0; i<HessianReduced.rows(); i++ )
        {
            HessianReduced(i,i) *= 1+Lambda;

        }

       // cout << "Hessian Red and Scaled:\n"<<Scale*HessianReduced<<endl<<endl;



        ErrorReduced = reducer.transpose()*Error;
        dU = reducer*Solver.compute(HessianReduced).solve(ErrorReduced);

//        if( std::isnan(dU.norm()) )
//        {
//            dU = reducer*Math::pseudoInverse(HessianReduced,ErrorReduced);
//        }

        //         Eigen::VectorXd uTmp(U.head(numControls+p_trajectory_definition->size()*constraint_list.size()));
        //         validate_jacobians( X_IC, uTmp, *p_trajectory_definition, constraint_list );
        //         if( constraint_list.size() && !constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N)) )
        //             int i=0;
        //         if( soft_constraint_list.size() && !soft_constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N)) )
        //             int i=0;


        U -= dU;

        // Enforce Control Constraints
        for(unsigned int j=0; j<p_trajectory_definition->size(); j++ )
        {
            double time  = (*p_trajectory_definition)[j].t_span.tf;
            Eigen::Matrix<double,M,1> state = (*p_trajectory_definition)[j].p_dynamic_system->x();

            for( unsigned int i=0; i < constraint_list.size(); i++ )
            {
                if(constraint_list[i]->has_u_sens() )
                    U.segment(j*N,N) = constraint_list[i]->inforce_u_constraint(time,state,U.segment(j*N,N));

            }
        }


        Eigen::VectorXd ErrorNew;
        costJacobian(U,ErrorNew);

        double deltaError = ErrorNew.norm() - Error.norm();

        bool initPass = true;
        while( deltaError > Error.norm()/100.0 && solveTimer.lap() < solveTime && Error.norm()/Error.size() <= nr_tol )
        {
            if(initPass)
            {
                Lambda *= 2.0;
                initPass = false;
            }

            //          cout << "Error Last: " << Error.norm() << " Error New: " << ErrorNew.norm() <<  " Delta Error: " << deltaError<<  " Lambda: " << Lambda << endl;


            double lambda_factor = (1.0+Lambda)/(1.0+Lambda/2.0);
            for( int i = 0; i<HessianReduced.rows(); i++ )
            {
                HessianReduced(i,i) *= lambda_factor;
                //Scale(i,i) = 1.0/HessianReduced.col(i).norm();
            }

            //cout << "Scale:\n" << HessianReduced <<endl<<endl;

            // dU = reducer*Scale*Solver.compute(HessianReducedSq*Scale).solve(ErrorReduced);
            dU = reducer*Solver.compute(HessianReduced).solve(ErrorReduced);

            U = U_last - dU;

            // cout << "dU: " << (dU).transpose() << endl <<endl;
            // cout << "Hess:\n" << Hessian << endl <<endl;
            // cout << "Hess RSq:\n"<< HessianReducedSq <<"\n------------"<< endl << endl;


            // Enforce Control Constraints
            for(unsigned int j=0; j<p_trajectory_definition->size(); j++ )
            {
                double time  = (*p_trajectory_definition)[j].t_span.tf;
                Eigen::Matrix<double,M,1> state = (*p_trajectory_definition)[j].p_dynamic_system->x();

                for( unsigned int i=0; i < constraint_list.size(); i++ )
                {
                    if(constraint_list[i]->has_u_sens() )
                        U.segment(j*N,N) = constraint_list[i]->inforce_u_constraint(time,state,U.segment(j*N,N));
                }
            }

            costJacobian(U,ErrorNew);
            deltaError = ErrorNew.norm() - Error.norm();


//            if( Lambda > 100 )
//            {
//                double solveTime = solveTimer.lap();

//                Eigen::VectorXd uTmp(U.head(numControls+p_trajectory_definition->size()*constraint_list.size()));
//                validate_jacobians( X_IC, uTmp, *p_trajectory_definition, constraint_list, soft_constraint_list,true );
//                if( constraint_list.size() && !constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N)) )
//                    int i=0;
//                if( soft_constraint_list.size() && !soft_constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N)) )
//                    int i=0;

//                cout << "Final: Error Last: " << Error.norm() << " Iterations: " << iterations <<  " Lambda: " << Lambda << " Solve Time: " << solveTime << endl<<endl;
//            }


            //            cout << "Error: " << Error.norm() <<" Delta Error: " << deltaError << " Lambda: " << Lambda << endl;
            Lambda *= 2.0;
        }
        if(deltaError < 0)
            Lambda /= 2.0;

        Lambda = std::max(Lambda,1e-10);

        //        cout << "Error Last: " << Error.norm() << " Error New: " << ErrorNew.norm() <<  " Delta Error: " << deltaError<< " Lambda: " << Lambda << " Solve Time: " << solveTimer.lap() << endl;
        Error = ErrorNew;
        U_last = U;



        //U -= Solver.compute(Hessian).solve(Error);
        iterations ++;

 //        cout << "Iteration: " << iterations << "\tError: "<< Error.norm() << endl;

    }while(Error.norm()/Error.size()>nr_tol && solveTimer.lap() < solveTime);



    //    if( solveTimer.lap() >= solveTime )
    //    {
    //        double solveTime = solveTimer.lap();

    //        Eigen::VectorXd uTmp(U.head(numControls+p_trajectory_definition->size()*constraint_list.size()));
    //        validate_jacobians( X_IC, uTmp, *p_trajectory_definition, constraint_list, soft_constraint_list,true );
    //        if( constraint_list.size() && !constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N)) )
    //            int i=0;
    //        if( soft_constraint_list.size() && !soft_constraint_list[0]->verify_constraint_jacobians(0,p_trajectory_definition->back().p_dynamic_system->x(),U.segment(numControls-N,N)) )
    //            int i=0;

    //        cout << "Final: Error Last: " << Error.norm() << " Iterations: " << iterations <<  " Lambda: " << Lambda << " Solve Time: " << solveTime << endl<<endl;
    //    }

 //   cout << "Final: Error Last: " << Error.norm() << " error size: " << Error.size() << " converge Value: " << Error.norm()/Error.size() << " converge tol: " << nr_tol << " Iterations: " << iterations <<  " Lambda: " << Lambda << " Solve Time: " << solveTimer.lap() << endl;

    // if( iterations == 150 || solveTimer.lap() >= solveTime )
    // {
    //     cout << "No Solution Found!" << endl;
    //     //        nr_solver.solve(U_last, nr_tol, 100,true,LEVENBERG_MARQUARDT,1);
    //     //        U = nr_solver.solution();
    // }
    //cout << "Iterations " << iterations << endl;

    return Error.norm()/Error.size() < nr_tol;
}

template<int M, int N>
bool MPC_Solver<M,N>::validate_jacobians( const Eigen::Matrix<double,M,1>& Xo_, Eigen::VectorXd& control_nominal, std::vector< Trajectory<M,N> >& trajectory_definition, const std::vector< Dynamic_System_Constraint<M,N>* >& constraint_list, const std::vector< Dynamic_System_Constraint<M,N>* >& soft_constraint_list_, bool printDiff )
{
    validation_mode = true;
    typename std::vector< Trajectory<M,N> >::iterator it; // trajectory iterator

    std::vector< Dynamic_System_Constraint<M,N>* > constraint_list_old = this->constraint_list;
    std::vector< Dynamic_System_Constraint<M,N>* > soft_constraint_list_old = this->soft_constraint_list;

    std::vector< Trajectory<M,N> >* p_trajectory_definition_old = this->p_trajectory_definition;
    Eigen::VectorXd control_nominal_old = this->U;
    Eigen::VectorXd U_Target_old = this->U_Target;
    Eigen::VectorXd IC_old = this->X_IC;



    this->constraint_list = constraint_list;
    soft_constraint_list = soft_constraint_list_;



    X_IC = Xo_;
    // save trajectory definition pointer for use in other internal functions
    p_trajectory_definition = &trajectory_definition;

    int tmp_numConstraints = numConstraints;
    int tmp_numControls = numControls;
    int tmp_numOutputs = numOutputs;
    int tmp_numStates = numStates;

    numConstraints  = trajectory_definition.size() * this->constraint_list.size();
    numControls = trajectory_definition.size() *N; // no control for the IC
    numStates = trajectory_definition.size() *M; // no control for the IC

    // Initialize System Control inputs to be target controls
    U = control_nominal;//.setZero(numControls,1);
    U_Target.setZero(numControls,1);
    Ju.setZero(numStates, numControls);
    Jx.setIdentity(numStates, numStates);


    nontrivial_Yu = false;
    nontrivial_Yxx = false;
    nontrivial_S = false;
    nontrivial_Yuu = false;
    nontrivial_Yux = false;
    nontrivial_H = false;
    nontrivial_Cux = false;

    numOutputs = 0;
    it = trajectory_definition.begin();
    for(int rowcountTarget=0; it!= p_trajectory_definition->end(); it++,rowcountTarget+=N)
    {
        if( it->includeTargetControlInCost )//for kalman filtering
        {
            U_Target.block(rowcountTarget,0,N,1) = it->target_control;
        }
        numOutputs += it->p_dynamic_system->getNumOutputs();
        nontrivial_Yu |= it->p_dynamic_system->has_out_u_sense();
        nontrivial_Yxx |= it->p_dynamic_system->has_out_xx_sense();
        nontrivial_S |= it->p_dynamic_system->has_out_xx_sense() || it->p_dynamic_system->has_out_ux_sense() || it->p_dynamic_system->has_ode_xu_sense() || it->p_dynamic_system->has_ode_uu_sense();
        nontrivial_Yuu |= it->p_dynamic_system->has_out_uu_sense();
        nontrivial_Yux |= it->p_dynamic_system->has_out_ux_sense();
        nontrivial_H |= it->p_dynamic_system->has_out_xx_sense() || it->p_dynamic_system->has_ode_xx_sense();

    }


    if( this->constraint_list.size() )
    {
        //typename std::vector< Dynamic_System_Constraint<M,N>* >::iterator it_c; // trajectory iterator

        for(unsigned int cons_ind = 0; cons_ind < this->constraint_list.size(); cons_ind ++ )//it_c = this->constraint_list.begin(); it_c!= this->constraint_list.end(); it_c++)
        {
            Dynamic_System_Constraint<M,N>* dsc = this->constraint_list[cons_ind];
            nontrivial_H |= dsc->has_xx_sens();
            nontrivial_Cux |= dsc->has_ux_sens();
        }
    }

    if( this->soft_constraint_list.size() )
    {
        //typename std::vector< Dynamic_System_Constraint<M,N>* >::iterator it_c; // trajectory iterator

        for(unsigned int cons_ind = 0; cons_ind < this->soft_constraint_list.size(); cons_ind ++ )//it_c = this->constraint_list.begin(); it_c!= this->constraint_list.end(); it_c++)
        {
            Dynamic_System_Constraint<M,N>* dsc = this->soft_constraint_list[cons_ind];
            nontrivial_H |= dsc->has_xx_sens();
            nontrivial_Cux |= dsc->has_ux_sens();
        }
    }

    if(nontrivial_Cux)
        Cux.setZero(numControls,numStates);

    //Set Q and R matrices to zeros
    Q.setZero(numOutputs,numOutputs);
    R.setZero(numControls,numControls);


    Jg.setZero(numOutputs, numStates);
    E.setZero(1,numOutputs);

    Z.setZero(1,numStates);            //Corrected from being numOutput to numState
    Y.setZero(1,numOutputs);

    //Packing Q
    it = p_trajectory_definition->begin();
    for(int diagElement=0; it!= p_trajectory_definition->end(); it++)
    {
        int stepSumOutput = it->p_dynamic_system->getNumOutputs();
        Q.block(diagElement,diagElement,stepSumOutput,stepSumOutput) = it->Q;
        diagElement+=it->p_dynamic_system->getNumOutputs();
    }

    //Packing R
    it = p_trajectory_definition->begin();
    for(int diagElement=0; it!= p_trajectory_definition->end(); it++,diagElement+=N)
    {
        R.block(diagElement,diagElement,N,N) = it->R;
    }

    H.setZero(numStates,numStates);
    S.setZero(numControls,numControls);
    Joo.setZero(numControls,numControls);
    Jox.setZero(numControls,numStates);


    // cout << U_Target.transpose() << endl;

    assert(control_nominal.rows() == numControls+numConstraints && "Miss Match in number of nominal controls, dont forget the lambda's for the constraints");
    //<- Verifies the jacobians
    double cost_nom = cost(control_nominal);
    Eigen::VectorXd costJacobianNominal,costJacobianNumeric,tmp;
    Eigen::MatrixXd costHessianNominal,costHessianNumeric;

    costJacobian(control_nominal,costJacobianNominal);
    costHessian(control_nominal,costHessianNominal);

    costHessianNumeric.setZero(numControls+numConstraints, numControls+numConstraints);
    costJacobianNumeric.setZero(numControls+numConstraints,1);

    Eigen::MatrixXd JxJu_nom = JxJu;
    Eigen::MatrixXd JxJu_num = JxJu;
    JxJu_num.setZero(numStates,numControls);

    Eigen::MatrixXd JgJxJu_nom = Jg*JxJu;
    Eigen::MatrixXd JgJxJu_num = Jg*JxJu;
    JgJxJu_num.setZero(numOutputs,numControls);

    Eigen::VectorXd x_nom;
    Eigen::VectorXd y_nom;
    x_nom.setZero(numStates,1);
    y_nom.setZero(numStates,1);
    int trajInd = 0;
    for( it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it ++ )
    {
        x_nom.block(trajInd*M,0,M,1) = it->p_dynamic_system->x();
        y_nom.block(trajInd*M,0,M,1) = it->p_dynamic_system->y();
        trajInd++;
    }



    for(int i=0; i<control_nominal.rows(); i++ )
    {
        double delta = std::max(abs(control_nominal(i))*1e-7,1e-7);
        control_nominal(i) += delta;
        costJacobianNumeric(i) = cost(control_nominal);

        costJacobian(control_nominal,tmp);
        costHessianNumeric.block(0,i,numConstraints+numControls,1) = tmp;

        int trajInd = 0;
        int outputNum = 0;
        if( i < numControls )
            for( it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it ++ )
            {
                JxJu_num.block(trajInd*M,i,M,1) = it->p_dynamic_system->x();
                JgJxJu_num.block(outputNum,i,it->p_dynamic_system->getNumOutputs(),1) = it->p_dynamic_system->getOutSelectMatrix()*it->p_dynamic_system->y();
                trajInd++;
                outputNum+=it->p_dynamic_system->getNumOutputs();
            }


        control_nominal(i) -= 2.0*delta;

        costJacobianNumeric(i) -= cost(control_nominal);
        costJacobianNumeric(i) /= 2.0*delta;

        costJacobian(control_nominal,tmp);
        costHessianNumeric.block(0,i,numConstraints+numControls,1) -= tmp;
        costHessianNumeric.block(0,i,numConstraints+numControls,1) /= 2.0*delta;


        trajInd = 0;
        outputNum = 0;
        if( i < numControls )
            for( it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it ++ )
            {
                JxJu_num.block(trajInd*M,i,M,1) -= it->p_dynamic_system->x();
                JgJxJu_num.block(outputNum,i,it->p_dynamic_system->getNumOutputs(),1) -= it->p_dynamic_system->getOutSelectMatrix()*it->p_dynamic_system->y();
                trajInd++;
                outputNum+=it->p_dynamic_system->getNumOutputs();
            }


        trajInd = 0;
        outputNum = 0;
        if( i < numControls )
            for( it = p_trajectory_definition->begin(); it!= p_trajectory_definition->end(); it ++ )
            {
                JxJu_num.block(trajInd*M,i,M,1) /= 2.0*delta;
                JgJxJu_num.block(outputNum,i,it->p_dynamic_system->getNumOutputs(),1) /= 2.0*delta;
                trajInd++;
                outputNum+=it->p_dynamic_system->getNumOutputs();
            }


        control_nominal(i) += delta;
    }

    double JxJu_error = (JxJu_nom - JxJu_num).norm()/JxJu_nom.norm();
    double JgJxJu_error = (JgJxJu_nom - JgJxJu_num).norm()/JgJxJu_nom.norm();
    double jacobian_error = (costJacobianNominal - costJacobianNumeric).norm()/costJacobianNominal.norm();
    double hessian_error = (costHessianNominal-costHessianNumeric).norm()/costHessianNominal.norm();

    if( JxJu_nom.norm() < 1e-3 )
        JxJu_error = (JxJu_nom - JxJu_num).norm();

    if( JgJxJu_nom.norm() < 1e-3 )
        JgJxJu_error = (JgJxJu_nom - JgJxJu_num).norm();

    if( costJacobianNominal.norm() < 1e-3 )
        jacobian_error = (costJacobianNominal - costJacobianNumeric).norm();

    if( costHessianNominal.norm() < 1e-3 )
        hessian_error = (costHessianNominal - costHessianNumeric).norm();


    cout << "Nominal Cost: " <<  cost_nom << endl;

    cout << "JxJu Error: " << JxJu_error*100 << "%" << endl;
    cout << "JgJxJu Error: " << JgJxJu_error*100 << "%" << endl;
    cout << "Jacobian Error: " << jacobian_error*100 << "%" << endl;
    cout << "Hessian Error: " <<  hessian_error*100 << "%" << endl;
    cout << "U: " << (control_nominal_old.head(numControls)+U_Target_old).transpose() << endl;
    double tol = 1e-7;

    if( printDiff && (JxJu_error > tol || isnan(JxJu_error)))
    {
        //cout << "Jx:\n" << Jx << endl << endl;
        //cout << "Ju:\n" << Ju << endl << endl;
        cout << "JxJu Numeric:\n" << JxJu_num << endl<< endl;
        cout << "JxJu Analitic:\n" << JxJu_nom << endl << endl;
        cout << "JxJu Difference:\n" << JxJu_nom - JxJu_num << endl<< endl;

    }

    if( printDiff && (JgJxJu_error > tol || isnan(JgJxJu_error)))
    {
        cout << "JgJxJu Numeric:\n" << JgJxJu_num << endl<< endl;
        cout << "JgJxJu Analitic:\n" << JgJxJu_nom << endl << endl;
        cout << "JgJxJu Difference:\n" << JgJxJu_nom - JgJxJu_num << endl<< endl;

    }

    if( printDiff && (jacobian_error > tol || isnan(jacobian_error)))
    {
        cout << "Jacobian Numeric: " << costJacobianNumeric.transpose() << endl;
        cout << "Jacobian Analitic: " << costJacobianNominal.transpose() << endl;
        cout << "Jacobian Difference: " << (costJacobianNominal - costJacobianNumeric).transpose() << endl;

        //        cout << "JxJu Numeric:\n" << JxJu_num << endl<< endl;
        //        cout << "JxJu Analitic:\n" << JxJu_nom << endl << endl;
        //        cout << "JxJu Difference:\n" << JxJu_nom - JxJu_num << endl<< endl;

        //        cout << "Hessian Numeric:\n" << costHessianNumeric << endl;
        //        cout << "Hessian Analitic:\n" << costHessianNominal << endl;
        //        cout << "Hessian Difference:\n" << (costHessianNominal - costHessianNumeric) << endl;

        cout << "------" << endl;
    }

    if( printDiff && (hessian_error > tol || isnan(hessian_error)))
    {
        cout << "Hessian Numeric:\n" << costHessianNumeric << endl;
        cout << "Hessian Analitic:\n" << costHessianNominal << endl;
        cout << "Hessian Difference:\n" << (costHessianNominal - costHessianNumeric) << endl;

    }


    this->soft_constraint_list = soft_constraint_list_old;
    this->constraint_list = constraint_list_old;
    this->p_trajectory_definition =  p_trajectory_definition_old;
    this->U = control_nominal_old;
    this->X_IC =  IC_old;
    this->U_Target = U_Target_old;

    numConstraints = tmp_numConstraints;
    numControls = tmp_numControls;
    numOutputs = tmp_numOutputs;
    numStates = tmp_numStates;

    validation_mode = false;

    return false;
}


#endif
