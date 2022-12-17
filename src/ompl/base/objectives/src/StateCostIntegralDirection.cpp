/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, Rice University
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Rice University nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/
  
 /* Author: Luis G. Torres */
  
 //#include "StateCostIntegralDirection.h"
 #include "ompl/base/spaces/RealVectorStateSpace.h"
 #include "ompl/base/objectives/StateCostIntegralDirection.h"
 #include <ctime>
 #include <fstream>


/* ompl::base::StateCostIntegralDirection::StateCostIntegralDirection(const SpaceInformationPtr &si,
                                                                    ompl::base::CostMatrix *cost_matrix,
                                                                    bool enableMotionCostInterpolation,
                                                                    float min_speed, float max_speed)
   : OptimizationObjective(si), interpolateMotionCost_(enableMotionCostInterpolation)
 {
     description_ = "State Cost Integral";     
     this->cost_matrix = cost_matrix;
     this->min_speed = min_speed;
     this->max_speed = max_speed;
 }
  */

ompl::base::StateCostIntegralDirection::StateCostIntegralDirection(const SpaceInformationPtr &si,
                                                                    const Eigen::Ref<Eigen::MatrixXd> &A,
                                                                    bool enableMotionCostInterpolation,
                                                                    float min_speed, float max_speed)
   : OptimizationObjective(si), interpolateMotionCost_(enableMotionCostInterpolation)
 {
     description_ = "State Cost Integral";     
     this->cost_matrix.loadMatrix(A);
     this->min_speed = min_speed;
     this->max_speed = max_speed;
 }
 ompl::base::Cost ompl::base::StateCostIntegralDirection::stateCost(const State *s) const
 {
    const double *pos = s->as<ompl::base::RealVectorStateSpace::StateType>()->values;
    double val = cost_matrix.interp(pos[0], pos[1]);
    if (isnan(val))
        return Cost(10e300);
    else{
         return Cost(val);
    }
 }
  
double ompl::base::StateCostIntegralDirection::ratioCost(const State *s1, const State *s2) const {
    const double *pos1 = s1->as<ompl::base::RealVectorStateSpace::StateType>()->values;
    const double *pos2 = s2->as<ompl::base::RealVectorStateSpace::StateType>()->values;

    double diff[2];
    for (int i = 0; i < 2; i++){
        diff[i] = pos2[i] - pos1[i];
    }
    double ratio = diff[0]/diff[1];

    if (ratio > min_speed && ratio < max_speed && diff[0] > 0 && diff[1] > 0){
        return 0.001;
    }    
    else if(diff[0] > 0 && diff[1] > 0){
        return pow(abs(ratio-1)+1, 4);
    }
    else{
        return 10e300;
        //return std::numeric_limits<float>::infinity();
    }

}

 ompl::base::Cost ompl::base::StateCostIntegralDirection::motionCost(const State *s1, const State *s2) const
 {
     if (interpolateMotionCost_)
     {
         Cost totalCost = this->identityCost();
  
         int nd = si_->getStateSpace()->validSegmentCount(s1, s2);

        const double *pos1 = s1->as<ompl::base::RealVectorStateSpace::StateType>()->values;
        const double *pos2 = s2->as<ompl::base::RealVectorStateSpace::StateType>()->values;
        double diff[2];
        for (int i = 0; i < 2; i++){
            diff[i] = pos2[i] - pos1[i];
        }        
        nd = round(sqrt(pow(diff[0], 2)+pow(diff[1], 2)));


        
         State *test1 = si_->cloneState(s1);
         Cost prevStateCost = this->stateCost(test1);

         if (nd > 1)
         {
             State *test2 = si_->allocState();
             for (int j = 1; j < nd; ++j)
             {                
                 si_->getStateSpace()->interpolate(s1, s2, (double)j / (double)nd, test2);
                 Cost nextStateCost = this->stateCost(test2);
                 totalCost = Cost(totalCost.value() +
                                  this->trapezoid(prevStateCost, nextStateCost, si_->distance(test1, test2)).value());
                 std::swap(test1, test2);
                 prevStateCost = nextStateCost;
             }
             si_->freeState(test2);
         }
  
         // Lastly, add s2
         double direction_cost = this->ratioCost(s1, s2);
//         std::cout << direction_cost << std::endl;
         totalCost = Cost((totalCost.value() +
                          this->trapezoid(prevStateCost, this->stateCost(s2), si_->distance(test1, s2)).value())*direction_cost);
  
         si_->freeState(test1);
  
         return totalCost;
     }
         double direction_cost = this->ratioCost(s1, s2);  
         return Cost(this->trapezoid(this->stateCost(s1), this->stateCost(s2), si_->distance(s1, s2)).value() + direction_cost*si_->distance(s1, s2));
 }
  
 bool ompl::base::StateCostIntegralDirection::isMotionCostInterpolationEnabled() const
 {
     return interpolateMotionCost_;
 }


void ompl::base::CostMatrix::loadCsv(const std::string &path){
    std::vector< std::vector<double>::iterator > grid_iter_list;
    std::vector<double> grid1, grid2, values;
    array<int,2> grid_sizes;    
    std::ifstream indata;
    indata.open(path);
    std::string line;
    int rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    for (int i = 0; i < rows; i++){
        grid1.push_back(i);
        grid2.push_back(i);
    }
    grid_iter_list.clear();
    grid_iter_list.push_back(grid1.begin());
    grid_iter_list.push_back(grid2.begin());    
    grid_sizes[0] = rows;
    grid_sizes[1] = rows;
    interp_ML = new InterpMultilinear<2, double>(grid_iter_list.begin(), grid_sizes.begin(), values.data(), values.data()+rows*rows);    
}

void ompl::base::CostMatrix::loadMatrix(const Eigen::Ref<Eigen::MatrixXd> &A){
    std::vector< std::vector<double>::iterator > grid_iter_list;
    std::vector<double> grid1, grid2, values;
    array<int,2> grid_sizes;    

//    grid_size = rows;
    for (int i = 0; i < A.rows(); i++){
        grid1.push_back(i);
    }
    for (int i = 0; i < A.cols(); i++){
        grid2.push_back(i);
    }
    grid_iter_list.clear();
    grid_iter_list.push_back(grid1.begin());
    grid_iter_list.push_back(grid2.begin());    
    grid_sizes[0] = A.rows();
    grid_sizes[1] = A.cols();
    interp_ML = new InterpMultilinear<2, double>(grid_iter_list.begin(), grid_sizes.begin(), A.data(), A.data()+A.rows()*A.cols());    

}

double ompl::base::CostMatrix::interp(const double x, const double y) const{
    array<double, 2> args = {x, y};
    return interp_ML->interp(args.begin());
}
