#ifndef MPC_SOLVER_H
#define MPC_SOLVER_H

#include <vector>
#include <Eigen/Dense>
//#include "Utilities/Math/NewtonRaphson.h"
#include "dynamic_system.h"
#include "dynamic_system_constraint.h"
#include "trajectory.h"
#include "tspan.h"
#include "r3tensor.h"
#include "Timer.h"

template<int M, int N>
class MPC_Solver
{
private:
  //  NewtonRaphson<MPC_Solver> nr_solver;


    Timing::stopWatch solveTimer;
    double solveTime;

    /* This is a helper function for the newton raphson class.
     *  It shoudl calculate the derivative of our cost function (what I'll call our error vector for the purposes of optimization)
     *  as a function of the current_controls.  Newton_raphson normally calls this function first. So we shoudl update (resolve)
     *  the system each time its called.
     */
    void costJacobian( Eigen::VectorXd& current_controls, Eigen::VectorXd& return_costJacobian);

    /* This is a helper function for the newton raphson class.
     *  It shoudl calculate the error-jacobian for the specified control values.
     *  NewtonRaphson always(normally) calls calculateError first, so we shoudl check to see
     *  if the control has changed before re-solving everything!  This (nominally) shoudl just be
     *  a packing function for us.
     */
    void costHessian( Eigen::VectorXd& current_controls, Eigen::MatrixXd& return_costHessian );

    /* This is a helper function to calculate the current cost for use in numerical verificaion
     * of the costJacobian and costHessian functions.
     */
    double cost(  const Eigen::VectorXd& current_controls );


    std::vector< Trajectory<M,N> >* p_trajectory_definition; ///<- a pointer to the trajectory defintion for use in the internal functions
    std::vector< Dynamic_System_Constraint<M,N>* > constraint_list; ///<- a pointer to the constraint list for use in the internal functions
    std::vector< Dynamic_System_Constraint<M,N>* > soft_constraint_list; ///<- a pointer to the soft constraint list for use in the internal functions

protected:
    //Eigen::MatrixXd C; // output selection matrix
    //Eigen::MatrixXd CtQC;
    Eigen::MatrixXd Jg;
    Eigen::MatrixXd Jo;
    Eigen::MatrixXd Joo;
    Eigen::MatrixXd Jox;
    Eigen::MatrixXd JgJxJu;
    Eigen::MatrixXd JoxJxJu;
    Eigen::MatrixXd Jx;
    Eigen::MatrixXd Ju;
    Eigen::MatrixXd JxJu;
    Eigen::MatrixXd H;
    Eigen::MatrixXd S;
    Eigen::MatrixXd S_c;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd R;
    Eigen::RowVectorXd Z; // E'QJgJx
    Eigen::RowVectorXd Y; // E'Q
    Eigen::RowVectorXd E; // output - output des
    Eigen::VectorXd U; // current control vector
    Eigen::VectorXd U_Target;
    Eigen::Matrix<double,M,1> X_IC; // initial condition
    Eigen::MatrixXd Cux;

    Eigen::MatrixXd x_der_constraints;
    Eigen::MatrixXd u_der_constraints;
    Eigen::RowVectorXd Z_constraint; // sum(Lambda*dc/dx)
    Eigen::RowVectorXd Z_soft_constraint; // sum(Lambda*dc/dx)



    int numConstraints;
    int numControls;
    int numOutputs;
    int numStates;

    bool nontrivial_Yu;
    bool nontrivial_Yxx;
    bool nontrivial_Yuu;
    bool nontrivial_Yux;
    bool nontrivial_Cux;

    bool nontrivial_S;
    bool nontrivial_H;

    bool validation_mode;

    bool newton_raphson();
    double int_tol;
    double nr_tol;

public:
    MPC_Solver(double step_integration_tol = 5e-8, double newton_rapson_tol = 1e-4, double solveTime = 0.030);
    bool solve( const Eigen::Matrix<double,M,1>& Xo_, std::vector< Trajectory<M,N> >& trajectory_definition, const std::vector< Dynamic_System_Constraint<M,N>* >& hard_constraint_list = std::vector< Dynamic_System_Constraint<M,N>* >(), const std::vector< Dynamic_System_Constraint<M,N>* >& soft_constraint_list = std::vector< Dynamic_System_Constraint<M,N>* >()); ///<- solves MPC problem and saves solution back into the trajectory provided
    bool validate_jacobians( const Eigen::Matrix<double,M,1>& Xo_, Eigen::VectorXd& control_nominal, std::vector< Trajectory<M,N> >& trajectory_definition, const std::vector< Dynamic_System_Constraint<M,N>* >& constraint_list = std::vector< Dynamic_System_Constraint<M,N>* >(), const std::vector< Dynamic_System_Constraint<M,N>* >& soft_constraint_list = std::vector< Dynamic_System_Constraint<M,N>* >(), bool printDiff = true); ///<- solves MPC problem and saves solution back into the trajectory provided

    //void initialize(const Eigen::Matrix<double,M,1>& Xo_, std::vector< Trajectory<M,N> >& trajectory_definition ); ///<- solves MPC problem and saves solution back into the trajectory provided
};

#include "mpc_solver.hpp"

#endif // MPC_SOLVER_H
