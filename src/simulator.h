#include <memory>
#include "model.h"
#include "controller.h"

class Simulator{
public:
    Simulator(){}
    Simulator(std::shared_ptr<Model> model0,std::shared_ptr<Controller> controller0);
    ~Simulator(){}

    void capture_2d();
    void captureAndTranslate_2d_withCargo();

    
private:
    std::shared_ptr<Model> model;
    std::shared_ptr<Controller> controller;
    int assignmentFrequency;
    int controlFrequency;
};