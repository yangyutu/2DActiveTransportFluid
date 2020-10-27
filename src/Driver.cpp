#include<iostream>
#include "simulator.h"
#include "controller.h"
#include "model.h"
#include "common.h"

Parameter parameter;

void readParameter();

int main(){

    readParameter();
    std::shared_ptr<Model> model(new Model());
    std::shared_ptr<Controller> controller(new Controller(model->getCurrState(), model->getTargets()));
    Simulator simulator(model, controller);
    
    for (int i = 0; i < parameter.nCycles; i++) {
        if (parameter.captureOnlyFlag == 1) {
            simulator.capture_2d();
        } else if (parameter.captureAndTransportFlag == 1) {
            simulator.captureAndTranslate_2d_withCargo();
        } 
    }
    return 0;
}

void readParameter(){
    std::string line;
    std::ifstream runfile;
    runfile.open("run.txt");
    getline(runfile, line);
    runfile >> parameter.N;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.radius;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.nCycles;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.captureStep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.captureOnlyFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.captureAndTransportFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.CollectiveMoveCycle >> parameter.collective_MoveStep >> parameter.collective_RestoreStep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dt;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.controlTimeInterval >> parameter.assignmentTimeInterval;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.diffu_t;    
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.diffu_r;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.maxVelocity;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.Bpp;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.Os_pressure;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.L_dep;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.cutoff;   
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.kappa;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.seed;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.trajOutputInterval;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.dynamicTargetFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.targetDiffuseRatio >> parameter.targetVelocityRatio >> parameter.cargoInteractingFlag;
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.targetCenter[0] >> parameter.targetCenter[1] >> parameter.targetCenter[2];
    getline(runfile, line);
    getline(runfile, line);
    runfile >> parameter.transporter_nb_thresh >> parameter.transporter_dist_thresh >> parameter.transport_angle_thresh;
    getline(runfile, line);
    getline(runfile, line);
    getline(runfile, parameter.iniConfig);
    getline(runfile, line);
    getline(runfile, parameter.filetag);
    getline(runfile, line);
    getline(runfile, parameter.targetConfig);
}

