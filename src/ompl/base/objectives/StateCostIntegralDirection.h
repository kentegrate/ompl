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
  
 #ifndef OMPL_BASE_OBJECTIVES_STATE_COST_OPTIMIZATION_OBJECTIVE_DIRECTION_
 #define OMPL_BASE_OBJECTIVES_STATE_COST_OPTIMIZATION_OBJECTIVE_DIRECTION_
  
 #include "ompl/base/OptimizationObjective.h"
 #include "ompl/base/objectives/linterp.h"
 #include <Eigen/Core>
 #include <Eigen/Dense>
 #include <vector>  

 namespace ompl
 {
     namespace base
     {
        class CostMatrix{
        public:
            CostMatrix(){
            }
            void loadCsv(const std::string &path);
            void loadMatrix(const Eigen::Ref<Eigen::MatrixXd> &A);
            double interp(const double x, const double y) const;

            //int grid_size;
        

            InterpMultilinear<2, double>* interp_ML;    
        };         

         class StateCostIntegralDirection : public OptimizationObjective
         {
         public:
//             StateCostIntegralDirection(const SpaceInformationPtr &si, CostMatrix *cost_matrix, bool enableMotionCostInterpolation = false, float min_speed=0.66, float max_speed=1.5);
             StateCostIntegralDirection(const SpaceInformationPtr &si, const Eigen::Ref<Eigen::MatrixXd> &A, bool enableMotionCostInterpolation = false, float min_speed=0.66, float max_speed=1.5);
  
             Cost stateCost(const State *s) const override;
  
             Cost motionCost(const State *s1, const State *s2) const override;
  
             bool isMotionCostInterpolationEnabled() const;
             CostMatrix cost_matrix;
             float min_speed;
             float max_speed;
 

         protected:
             bool interpolateMotionCost_;
 
             Cost trapezoid(Cost c1, Cost c2, double dist) const
             {
                 return Cost(0.5 * dist * (c1.value() + c2.value()));
             }

             double ratioCost(const State *s1, const State *s2) const;
         };
     }
 }




  
 #endif