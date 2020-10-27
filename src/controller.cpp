#include<cmath>
#include<set>
#include<limits>
#include "controller.h"
#include "HungarianAlg.h"
#include<fstream>
#include "common.h"
#define M_PI 3.1415926
extern Parameter parameter;

Controller::Controller(Model::state s,Model::state targets){
    s_ = s;
    targets_ = targets;
    radius = parameter.radius;
    numP = parameter.N;
    numTargets = targets.size();
    dimP = 2;
}


void Controller::calControl2d(Model::state s, Model::state targets){
// here calculate the control option given the targets
    for(int i=0; i < numP; i++){
        s[i]->u = 0;

        int t_idx = s[i]->targetIdx;
        if (t_idx >= 0){ // t_idx >=0 indicates a valid target
            double rx, ry;
            double pos_x, pos_y; // particle position
            pos_x = (s[i]->r[0]) / radius;
            pos_y = (s[i]->r[1]) / radius;

            if (parameter.dynamicTargetFlag) {
                // if the target is dynamic, then we track its next mean position
                rx = s[i]->targetPos[0] + parameter.maxVelocity * parameter.targetVelocityRatio / radius
                    - pos_x;
                ry = s[i]->targetPos[1] - pos_y;
            }
            else {
                rx = s[i]->targetPos[0] - pos_x;
                ry = s[i]->targetPos[1] - pos_y;
            }

            double dot_prod = cos(s[i]->phi) * rx + sin(s[i]->phi) * ry;

            if (dot_prod < 0) {
                s[i]->u = 0;
            } else {
                
               double disp_limit =  parameter.controlTimeInterval*parameter.maxVelocity/radius;

                if (dot_prod > disp_limit) {
                    s[i]->u = 1;                    
                } else {
                    s[i]->u = dot_prod/disp_limit;
                }              
 
            }

        } 
    }
}



double Controller::calAssignment(Model::state s, Model::state targets) {
    
    // cost is a numP by numP matrix, recording all the pair Euclidean distances
    vector< vector<double> > Cost(numP, vector<double>(numP));
    double totalCost = 0.0;
    for(int i=0; i<numP; i++){
	for(int j=0; j<numP; j++){
            double rx = targets[j]->r[0] - s[i]->r[0]/radius ;
            double ry = targets[j]->r[1] - s[i]->r[1]/radius ;
            double rz = targets[j]->r[2] - s[i]->r[2]/radius ;
            double dist = rx*rx + ry*ry + rz*rz;
            dist = sqrt(dist);
            
            Cost[i][j] = dist;

        }
    }   
    AssignmentProblemSolver APS;
    APS.Solve(Cost, assignment);
    for(int i=0; i < numP; i++){
        s[i]->targetIdx = assignment[i];
        s[i]->targetPos[0] = targets[s[i]->targetIdx]->r[0];
        s[i]->targetPos[1] = targets[s[i]->targetIdx]->r[1];
        s[i]->targetPos[2] = targets[s[i]->targetIdx]->r[2];
        s[i]->cost = Cost[i][assignment[i]];
        s[i]->EudDistToTarget = Cost[i][assignment[i]];
        targets[s[i]->targetIdx]->targetIdx = i;
        totalCost += s[i]->cost;
    }
    return totalCost;
}


void Controller::calEudDist(Model::state s){
    double dist;
    double r[3];
    int size = s.size();
    for (int i = 0; i < size; i++) {
            dist = 0.0;
            for (int k = 0; k < dimP; k++) {
                r[k] = s[i]->targetPos[k] - s[i]->r[k] / radius;
                dist += pow(r[k], 2.0);
            }
        dist = sqrt(dist);
        s[i]->EudDistToTarget;
    }
}

void Controller::translate_2d(double phi, Model::state s) {
    // first select a subset of particles that is in the right direction
    //  let them go until shape change too much
    this->calInlier(s);

    for (int i = 0; i < numP; i++) {
        double proj;
        proj = cos(s[i]->phi - phi);



        s[i]->transporterFlag = 0;
        //s[i]->u = 0.0;
        if (s[i]->nbcount >= parameter.transporter_nb_thresh &&
            s[i]->EudDistToTarget < parameter.transporter_dist_thresh) {

            if (proj > cos(parameter.transport_angle_thresh / 180.0 * M_PI)) {
                s[i]->transporterFlag = 1;
                s[i]->u = 1.0;
            }
        }
    }
}

void Controller::calInlier(Model::state s) {

    double r[3], dist;
    int size = s.size();
    for (int i = 0; i < size; i++) {
        s[i]->nbcount = 0;
        s[i]->inlier = 0;
    }
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            dist = 0.0;
            for (int k = 0; k < dimP; k++) {

                r[k] = (s[j]->r[k] - s[i]->r[k]) / radius;
                dist += pow(r[k], 2.0);
            }
            dist = sqrt(dist);
            if (dist < 2.5) {
                s[i]->nbcount++;
                s[j]->nbcount++;
            }

        }
    }
    for (int i = 0; i < size; i++) {
        if (s[i]->nbcount >= 3) {
            s[i]->inlier = 1;
        }
    }

}