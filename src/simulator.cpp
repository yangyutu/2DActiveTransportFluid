#include "simulator.h"
#include "common.h"

extern Parameter parameter;

Simulator::Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0):
                    model(model0),controller(controller0){
    controlFrequency = parameter.controlTimeInterval / model->dt();
    assignmentFrequency = parameter.assignmentTimeInterval/parameter.controlTimeInterval;
}

void Simulator::capture_2d(){
    double totalCost;
    model->createInitialState();
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets());
    controller->calControl2d(model->getCurrState(),model->getTargets());
    std::cout << totalCost << std::endl;
    for(int s=0; s < parameter.captureStep; s++){
        if ((s + 1) % assignmentFrequency == 0){
            totalCost = controller->calAssignment(model->getCurrState(),model->getTargets());
            std::cout << "cargo catch in shape forming " << s << "\t" <<totalCost << std::endl;
        }
        controller->calControl2d(model->getCurrState(),model->getTargets());        
        model->run(controlFrequency);
    }
}


void Simulator::captureAndTranslate_2d_withCargo(){
    double totalCost;   
    model->createInitialState();

    // first do a catch
    totalCost = controller->calAssignment(model->getCurrState(),model->getTargets());
    controller->calControl2d(model->getCurrState(),model->getTargets());
    std::cout << totalCost << std::endl;
    for(int s=0; s < parameter.captureStep; s++){
        if ((s+1)%assignmentFrequency == 0){
            totalCost = controller->calAssignment(model->getCurrState(),model->getTargets());
            std::cout << "catch step "<<s << "\t" <<totalCost << std::endl;
        }
        controller->calControl2d(model->getCurrState(),model->getTargets());        
        model->run(controlFrequency);
    }

    for(int c = 0; c < parameter.CollectiveMoveCycle; c++){
        for (int s = 0; s < parameter.collective_MoveStep; s++) {
            if ((s + 1) % assignmentFrequency == 0){
                totalCost = controller->calAssignment(model->getCurrState(),model->getTargets());
                std::cout << "motion cycle: "<< c << "\t" <<totalCost << std::endl;
            }
                // first all the particles are choosed to following assembly path
            controller->calControl2d(model->getCurrState(),model->getTargets());
            // part of particles are selected as transporters
            controller->translate_2d(0.0, model->getCurrState());
            model->run(controlFrequency);
        }
        
        for (int s = 0; s < parameter.collective_RestoreStep; s++) {
            if ((s + 1) % assignmentFrequency == 0) {
                totalCost = controller->calAssignment(model->getCurrState(), model->getTargets());
                std::cout << "recover step: "<< s << "\t"<<totalCost <<totalCost << std::endl;
            }
            controller->calControl2d(model->getCurrState(), model->getTargets());
            model->run(controlFrequency);
        }
    
    }
}

