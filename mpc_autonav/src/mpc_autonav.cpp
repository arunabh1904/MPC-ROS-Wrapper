#include <ros/ros.h>
#include <eigen3/Eigen/Core>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <assert.h>
#include <stdio.h>
#include "math.h"

//************************************* OpenCV includes***************

#include <opencv2/core.hpp>
#include <cv.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

//************************************* Message_types includes***************

#include <geometry_msgs/TransformStamped.h> //for odom to base_link transform
#include <geometry_msgs/Twist.h> //twist type for velocity
#include <nav_msgs/Odometry.h> //husky pose twist data from odom filtered
#include <velodyne_msgs/VelodyneScan.h> //lidar scan as lidar packet
#include <sensor_msgs/PointCloud2.h> //pcl2 format to extract xyz points
#include <std_msgs/Float64.h> //velocity command type
//#include <sensor_msgs/point_cloud2_iterator.h> //iterator for reading data
#include <std_msgs/Float32.h>
//#include <velodyne_pointcloud/point_types.h>//POINTXYZIR
#include <sensor_msgs/PointField.h>

//************************************* tf2 includes***************

#include <tf2/transform_datatypes.h>//quaternion functions
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include "tf2_geometry_msgs/tf2_geometry_msgs.h"

//************************************* tf includes***************

#include <tf/transform_listener.h>
#include <tf/tf.h>
#include <tf/LinearMath/Quaternion.h>
#include <tf/LinearMath/Matrix3x3.h>

//************************************* PCL includes***************
#include "pcl_conversions/pcl_conversions.h"
#include "pcl_ros/point_cloud.h"
#include "pcl/point_types.h"

//************************************* MPC includes***************

#include "control_saturation_constraint.h"
#include "husky_collision_constraint.h"
//#include "husky_potential_field.h"
#include "husky_system.h"
#include "mpc_solver.h"
#include "trajectory.h"

#define LIDAR_DISPLACEMENT_DISTANCE -0.125 //0.125m
#define OBSTACLE_FROM_GROUND 0.02 // 2 cm
#define WHEEL_RADIUS 0.17775

#define PRINT_CONVERGED
#define PRINT_DEBUG
//#define SHOW_OBS_IMAGE
//#define DONT_MOVE
//#define SOLVED_VELOCITY

class Callback
{
private:
  //struct for pose data
  struct Pose{
    double x;
    double y;
    double heading;
  };

  //struct for storing pointcloud2 data
  struct Pcl2Data{
    float x;
    float y;
    float z;
    float buffer;
    float intensity;
    short ring;
    char buffer2[10];
  };


  //initialize the listener
  tf::TransformListener listener;//initialization should never be in a callback
  //tf2_ros::Buffer tfBuffer;
  //tf2_ros::TransformListener listener(tfBuffer);


public:
  Eigen::Matrix<double, 3, 4> robot2World_Transform;//identity at first, keep updating for every frame as we go.

  tf::StampedTransform odom_base_transform;

  //robot2World_Transform.setIdentity();
  std::vector<Eigen::Vector2d> obstacle_list;//list for storing obstacle points

  std::vector<Eigen::Vector3d> world_obstacle_list;//obstacle points

  Pose pose;

  //to get current pose data
  bool acquire_robotpose();

  //allocate image objects
  cv::Mat filter_image;// = cv::Mat(1875, 16, CV_8UC1, cv::Scalar(0));// = filter_image(1875, 16, CV_8UC1, 0);
  cv::Mat index_image;
  cv::Mat obstacle_image;

  //callback for pose message
  //void pose_message_received_cb(const nav_msgs::Odometry::ConstPtr& pose_msg);

  //callback for lidar points
  void points_received_cb(const sensor_msgs::PointCloud2::ConstPtr& pcloud);

};



int main(int argc, char **argv)
{

  /*****   SETUP ROS STUFF *******/
  //initialize ros node
  ros::init(argc, argv, "mpc_autonav");

  //create node handle
  ros::NodeHandle node_handle;

  Callback call_back;
  //ROS_INFO("MAIN");

  // loop update rate
  double node_update_rate_Hz = 20;
  ros::Rate rate( node_update_rate_Hz );//setup ros rate (command line check: rostopic hz)

  //create velocity publisher object
  int velocity_command_que_size = 1;
  ros::Publisher velocity_cmdPublisher = node_handle.advertise<geometry_msgs::Twist>("husky_velocity_controller/cmd_vel", velocity_command_que_size);

  //create subscriber object for pose
  //ros::Subscriber pose_Subscriber = node_handle.subscribe("odometry/filtered", pose_sub_rate_hz,&  Callback::pose_message_received_cb, &call_back);//odometry/filtered husky topic which publishes odometry data as nav_msgs/Odometry

  // create subscriber object for LIDAR
  int lidar_buffer_que_size = 1;
  //initialize lidar scan subscriber
  ros::Subscriber points_Subscriber = node_handle.subscribe("/velodyne_points", lidar_buffer_que_size,&  Callback::points_received_cb, &call_back);//velodyne_points topic which publishes lidar data as sensor_msgs/PointCloud2

  //*****************************************************************




  /*****************************  Initialization of MPC Objects   *****************************/

  std::vector< Trajectory<3,2> > pathTrajectory;//trajectory object

  MPC_Solver<3,2> mpcSolver(1e-8, 5e-4, 1.0/node_update_rate_Hz*0.90);//mpc solver object - int_tol, nr_tol, sol_time = 0.045

  // create constraints
  Control_Saturation_Constraint<3,2> linVel_Constraint(0,1.0,-1.0);//lin vel saturation
  Control_Saturation_Constraint<3,2> angVel_Constraint(1, M_PI/6.0, -M_PI/6.0);//ang vel saturation
  //+25% husky initialization
  Husky_Collision_Constraint obstacles(147, 1, 1.25, 0.875, 1.0625, Eigen::Vector3d(0,0,0));//collision constraint object
  //Husky_Collision_Constraint obstacles(5.0, 1, 2.5, 2, 1.0625, Eigen::Vector3d(0,0,0));//collision constraint object
  // create constraint list
  std::vector< Dynamic_System_Constraint<3,2>* > constraintList;//hard constraint list
  std::vector< Dynamic_System_Constraint<3,2>* > soft_ConstraintList;//soft constraint list

  // append constraints to list
  //constraintList.push_back(& linVel_Constraint);
  //constraintList.push_back(& angVel_Constraint);
  constraintList.push_back(& obstacles);
  soft_ConstraintList.push_back(& obstacles);


  // initialize MPC Weighting
  Eigen::Matrix2d R;//control penalty matrix
  R << 1.0/std::pow(0.07,2), 0,
      0, 1.0/std::pow(0.05,2);

  Eigen::Matrix<double,1,1> Q;//error penalty matrix
  //Q.setIdentity(3,3);
  Q(0,0) = 1.0/std::pow(M_PI/16.0,2);

  Eigen::Matrix<double,1,3> C;//selection matrix
  C << 0,0,1; // pick off theta only

  // Create trajectory list
  int numPoints = 10;//num of traj points to solve in the future
  double trajectory_time_delta = 0.5;//(1.0/node_update_rate_Hz)*16.0; // 0.8 seconds

  double constant_target_velocity = 1;
#ifdef DONT_MOVE
  constant_target_velocity = 0;
#endif

  double testAngle = 0;

  /*pathTrajectory.push_back( Trajectory<3,2>(HuskySystem()) );
    pathTrajectory.back().target_output = Eigen::Matrix<double,1,1>(testAngle);//set to target ouput THIS NEEDS TO BE UPDATED IN LOOP
    pathTrajectory.back().target_control = Eigen::Vector2d(constant_target_velocity, 0); // set to target vel  THIS NEEDS TO BE UPDATED IN LOOP
    pathTrajectory.back().Q = Q/4.0;
    pathTrajectory.back().R = R/4.0;
    pathTrajectory.back().includeTargetControlInCost = false;
    pathTrajectory.back().t_span = T_Span(0,trajectory_time_delta/4.0);//0.05 secs
    pathTrajectory.back().p_dynamic_system->setOutSelectMatrix(C);
  */
  for( int i=0; i<numPoints; i++)
  {
    pathTrajectory.push_back( Trajectory<3,2>(HuskySystem()) );
    pathTrajectory.back().target_output = Eigen::Matrix<double,1,1>(testAngle);
    pathTrajectory.back().target_control = Eigen::Vector2d(constant_target_velocity,0);
    pathTrajectory.back().Q = Q;
    pathTrajectory.back().R = R;
    pathTrajectory.back().includeTargetControlInCost = false;
    pathTrajectory.back().t_span = T_Span(i*trajectory_time_delta,(i+1)*trajectory_time_delta);
    pathTrajectory.back().p_dynamic_system->setOutSelectMatrix(C);
  }

  //cout << "Target Control: " << pathTrajectory.front().target_control.transpose() << endl;

  //    //get lidar points and set new obstacles in world frame into constraint list.
  //    call_back.acquire_robotpose();

  //    call_back.create_obstacle_list(call_back.obstacle_list, call_back.robot2World_Transform);

  //    obstacles.updateStateAndObstacle(Eigen::Vector3d(0,0,0), call_back.obstacle_list);


  //write a service client thing for this.
  //goal(std::Nan, std::Nan, std::Nan, std::Nan, std::Nan);

  /*****************************  END Initialization of MPC Objects   *****************************/



  /*******************   BEGIN LOOP *********************/

  while(ros::ok())
  {

    // Get robot pose
    call_back.acquire_robotpose();

    // extract state from callback object for clarity
    Eigen::Vector3d currentState( call_back.pose.x, call_back.pose.y, call_back.pose.heading );

    // PRINT FOR DEBUG
    //cout<<"current state: "<<currentState<<endl;
    //for(int trajind=0;trajind<numPoints;trajind++)
    //{
    //cout<<"pathTrajectory.target_output:  "<<pathTrajectory[trajind].target_output<<endl;
    //cout<<"pathTrajectory.target_control:  "<<pathTrajectory[trajind].target_control<<endl;
    //}
    //cout<<"constraintList: "<<&constraintList<<endl;
    //cout<<"soft_ConstraintList: "<<&soft_ConstraintList<<endl;


    //call_back.create_world_obstacle_list();
#ifdef PRINT_DEBUG
    //cout << "world obs size: " << call_back.world_obstacle_list.size() << endl;
#endif
    obstacles.updateStateAndObstacle(currentState, call_back.world_obstacle_list);

#ifdef SHOW_OBS_IMAGE

    int obsDist = 50;
    cv::Mat obsIm(obsDist,obsDist,CV_8UC1,cv::Scalar(0));
    for( int obs_ind = 0; obs_ind < call_back.world_obstacle_list.size(); obs_ind ++ )
    {
      int xloc = call_back.world_obstacle_list[obs_ind].x()+obsDist/2;
      int yloc = call_back.world_obstacle_list[obs_ind].y()+obsDist/2;

#ifdef PRINT_DEBUG
      //cout << "obs loc: " << call_back.world_obstacle_list[obs_ind].transpose() << endl;
#endif
      if( xloc >=0 && xloc < obsDist && yloc >=0 && yloc < obsDist)
      {
        obsIm.at<uchar>(yloc,xloc) = 255;
      }

    }

    cv::resize(obsIm,obsIm,cv::Size(800,800));
    float scale = 800.0/obsDist;

    cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(scale*(call_back.pose.x+obsDist/2),scale*(call_back.pose.y+obsDist/2)), cv::Size2f(2*scale,scale), call_back.pose.heading*180.0/M_PI);
    cv::Point2f vertices[4];
    rRect.points(vertices);
    for (int i = 0; i < 4; i++)
      line(obsIm, vertices[i], vertices[(i+1)%4],cv::Scalar(255),1);

    cv::imshow("obs layout", obsIm);
#endif


    // GET ROBOT DESIRED STATE AND CONTROL

    // UPDATE TRAJECTORY DEFINITION

    //double running_velocity = 0.5;

    //pathTrajectory.back().target_control = Eigen::Vector2d(running_velocity, 0);



    //solve mpc
    bool converged = mpcSolver.solve( currentState,  pathTrajectory,  constraintList,  soft_ConstraintList);
#ifdef PRINT_CONVERGED
    if (!converged){
      cout << "converged: "<< converged << endl;
    }
#endif
    //  HARD ENFORCE VELOCITY CONSTRAINTS JUST IN CASE
    linVel_Constraint.inforce_u_constraint(0, Eigen::Matrix<double,3,1>::Zero() ,  pathTrajectory[0].soln_control);
    angVel_Constraint.inforce_u_constraint(0, Eigen::Matrix<double,3,1>::Zero() ,  pathTrajectory[0].soln_control);


    // SETUP CONTROL COMMANDS
    double linear_velocity = pathTrajectory[0].soln_control(0);//assigns vel in x dir
    double angular_velocity = pathTrajectory[0].soln_control(1);//assigns ang vel about z dir

#ifdef SOLVED_VELOCITY
    cout<<"soln linear vel: "<<linear_velocity<<endl;
    cout<<"soln angular vel: "<<angular_velocity<<endl;
#endif
    if( !converged )
    {
      //check if collision occurs
      bool colliding  = (0 != obstacles.constraint_cost(0,  pathTrajectory[0].soln_state,  pathTrajectory[0].soln_control, HARD));

      linear_velocity = 0;
      angular_velocity = 0;

      if( colliding )
      {
        // Zero Velocities
        linear_velocity = 0;
        angular_velocity = 0;
        // Zero Controls
        for(int trajInd = 0; trajInd< numPoints; trajInd++)
        {
          //set solved controls to zero
          pathTrajectory[trajInd].soln_control = Eigen::Vector2d(0,0);
        }
      }
    }

    // COPY CONTROLS BACKWARDS DOWN TRAJECTORY TO HELP WITH SOLVING NEXT TIME STEP
    for( int trajIndex = 1; trajIndex < numPoints; trajIndex++ )
    {
      pathTrajectory[trajIndex-1].soln_control =  pathTrajectory[trajIndex].soln_control;
    }

    //   cout<<"linear vel: "<<linear_velocity<<endl;
    //  cout<<"angular vel: "<<angular_velocity<<endl;

    //-----------------------------------------------------------------------------------------------


    //Define Velocity Command
    geometry_msgs::Twist velocity_command;//twist message to send velocity command.

    //set velocities to be sent to the husky
    velocity_command.linear.x =  linear_velocity;
    velocity_command.angular.z =  angular_velocity;

    //publish velocity command to husky.
    velocity_cmdPublisher.publish(velocity_command);

    // Call callback functions for the velocity command published data
    ros::spinOnce();

    //wait for next iteration
    rate.sleep();
  }//ros::ok
  //shuts down the node gracefully (ctrl-c)

  geometry_msgs::Twist velocity_command;//twist message to send velocity command.

  //set velocities to be sent to the husky
  velocity_command.linear.x =  0;
  velocity_command.angular.z =  0;

  //publish velocity command to husky.
  velocity_cmdPublisher.publish(velocity_command);

  // Call callback functions for the velocity command published data
  ros::spinOnce();



  ros::shutdown();
}//main


void Callback::points_received_cb(const sensor_msgs::PointCloud2::ConstPtr& pCloud)
{
  float distThreshold = 30; // in meters...
  int pointsPerRing = 500;

  //image object initialization
  filter_image = cv::Mat(16,pointsPerRing, CV_8UC1, cv::Scalar(255));
  //blurred_image = cv::Mat(pointsPerRing, 16, CV_8UC1, cv::Scalar(255));
  index_image = cv::Mat( 16, pointsPerRing, CV_16UC1, cv::Scalar(pointsPerRing*16));
  obstacle_image = cv::Mat(16,pointsPerRing, CV_8UC1, cv::Scalar(0));
  //clear obstacle list
  obstacle_list.clear();

  //cout<<"pcloud height: "<<pCloud->height<<endl;//1
  //cout<<"pcloud pointstep: "<<pCloud->point_step<<endl;//32
  //cout<<"pcloud frame id: "<<pCloud->header.frame_id;//frame - velodyne
  //cout<<"pcloud rowstep: "<<pCloud->row_step<<endl;//32*numPoints
  //cout<<"pcloud is_bigendian: "<<pCloud->is_bigendian;
  //cout<<"pcloud width: "<<pCloud->width;//numPoints

  for (int numPoints = 0; numPoints<pCloud->row_step; numPoints+=32)
  {
    //declare pcl2 pointer for data copy
    const Pcl2Data *pcl2Data = (Pcl2Data*)(&pCloud->data[numPoints]);
    //calculate r
    float radial_distance = std::sqrt(pcl2Data->x * pcl2Data->x + pcl2Data->y * pcl2Data->y + pcl2Data->z * pcl2Data->z);

    //cout << "Ring no: " << pcl2Data->ring << "\tz: " << pcl2Data->z << endl;

    float angle = atan2(pcl2Data->y, pcl2Data->x);

    //Azimuthal filtering + removing dummy points on lidar + restricting field of view to forward 180 degrees
    if(pcl2Data->z < (-0.45) || radial_distance < 0.001 || angle>M_PI/2.0 || angle<-M_PI/2.0)
    {
      //cout << "Ring no: " << pcl2Data->ring << "\tz: " << pcl2Data->z << endl;
      //cout << "\tx: " << pcl2Data->x <<"\ty: " << pcl2Data->y <<"\tz: " << pcl2Data->z << endl;
      continue;
    }

    //scale r to 0 to 255 for storing as image pixels
    //if( radial_distance > 30 )
    //int r = radial_distance/255.0*distThreshold;
    int r = radial_distance*255.0/distThreshold;
    if( r > 255 )
      r = 255;
    //else
    //    r = 30;


    int th = (M_PI+angle)*pointsPerRing/(2.0*M_PI);

    //cout << "Radial Distance: " << radial_distance << "\tRad_int: " << r << "\t Ring no: " << pcl2Data->ring <<  "\tTheta Ind: " << th << endl ;

    //filter_image.data[(filter_image.cols*pcl2Data->ring+(numPoints/32))] = r;
    filter_image.data[th+pointsPerRing*pcl2Data->ring]  = r;
    index_image.at<unsigned short>(pcl2Data->ring,th) = numPoints;

  }

  //Gaussian Blur - img size - 16*500
  //cv::GaussianBlur(filter_image,filter_image,cv::Size(5,5),3,3);

  world_obstacle_list.clear();
  int ignored = 0;
  for( int r = 1; r < filter_image.rows-1; r++ )
  {
    for ( int c = 1; c< filter_image.cols-1; c++ )
    {
      int index = index_image.at<unsigned short>(r,c);
      const Pcl2Data *pcl2Data = (Pcl2Data*)(&pCloud->data[index]);

      unsigned char thisVal = filter_image.data[c+filter_image.step*r];

      if( thisVal >= 250 )
      {
        // skip these points
        continue;
      }


      bool isMin = true;

      //for( int ud = -1; ud < 2 && isMin; ud += 2 )
      for( int ud = -1; ud < 2 && isMin; ud += 2 )// stop if its not the min and move on
      {
        for( int lr = -1; lr < 2 && isMin; lr += 2 )
        {
          isMin &= (thisVal < filter_image.data[(c+ud)+filter_image.step*(r+lr)]);
        }
      }

      if( isMin )
      {
        // push back point
        Eigen::Vector4d pclPoint = Eigen::Vector4d(pcl2Data->x,pcl2Data->y,pcl2Data->z,1);
        Eigen::Vector3d worldPoint = robot2World_Transform * Eigen::Vector4d(pcl2Data->x,pcl2Data->y,pcl2Data->z,1);
#ifdef PRINT_DEBUG
        //cout << "pcl point " << pclPoint.transpose() << "\tworld Point: " << worldPoint.transpose() << "\twp norm: " <<worldPoint.norm() << endl;
#endif
        world_obstacle_list.push_back(worldPoint);
        obstacle_image.data[(c)+obstacle_image.step*(r)] = 255;
        //obstacle_image.at<unsigned short>(r,c) = 255;
      } else
      {
        ignored ++;
      }
      //cout << "r: " << r << "\tc: " << c << endl;
    }
  }

#ifdef SHOW_OBS_IMAGE
  cv::imshow("filter_image", filter_image);
  cv::imshow("obstacles",obstacle_image);
  cv::waitKey(1);
#endif

#ifdef PRINT_DEBUG
  cout << "Number of obstacles: " << world_obstacle_list.size() << "\tNumber Ignored: " << ignored << endl;
#endif
}

bool Callback::acquire_robotpose()
{
  //geometry_msgs::TransformStamped odom_base_transform;//odom2basetransform
  //tf::StampedTransform odom_base_transform;
  tf::Matrix3x3 mat;
  double roll;
  double pitch;
  double yaw;

  try
  {
    // this is needed for first frame incase the EKF has not started yet
    listener.waitForTransform("odom", "base_link", ros::Time(0), ros::Duration(1.0));
    listener.lookupTransform("odom", "base_link", ros::Time(0), odom_base_transform);
  }
  catch(tf::TransformException ex)
  {
    ROS_ERROR("%s",ex.what());
    ROS_INFO("odom transform not received");
    return false;
  }

  mat.setRotation(odom_base_transform.getRotation());

  for( int i=0;i<3;i++ )
    for( int j=0;j<3;j++ )
      robot2World_Transform(i,j) = mat[i][j]; // copy rotation matrix over to eigen version

  // copy over vector... Maybe we should have just used tf types... :)
  robot2World_Transform.block(0,3,3,1) = Eigen::Vector3d(odom_base_transform.getOrigin().x(),odom_base_transform.getOrigin().y(),odom_base_transform.getOrigin().z()+WHEEL_RADIUS);


  //mat.getRPY(roll, pitch, yaw);

  pose.x = odom_base_transform.getOrigin().getX();
  pose.y = odom_base_transform.getOrigin().getY();
  pose.heading = atan2(robot2World_Transform(1,0),robot2World_Transform(0,0));

  //    double sth = std::sin(pose.heading);
  //    double cth = std::cos(pose.heading);

  //    robot2World_Transform << cth, -sth,  pose.x,
  //            sth,  cth,  pose.y;

  //cout<< "robot2World_Transform : "<<robot2World_Transform;
  //cout<<odom_base_transform;
  //cout<<"x: "<<pose.x<<endl;
  //cout<<"y: "<<pose.y<<endl;
  //cout<<"heading: "<<pose.heading<<endl;
  /*try
    {
        //velodyne to base_link
        listener.waitForTransform("base_link", "velodyne", ros::Time(0), ros::Duration(1.0));
        listener.lookupTransform("base_link", "velodyne", ros::Time(0), base_velo_transform);
    }
    catch(tf::TransformException ex)
    {
        ROS_ERROR("%s",ex.what());
        ROS_INFO("velo transform not received");
        return false;
    }

    velodyne_z = base_velo_transform.getOrigin().getZ();

    cout<<"velo_z: "<< velodyne_z << endl;*/

  return true;
}
