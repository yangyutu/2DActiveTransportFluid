#pragma once
#include<vector>
#include<map>
#include <random>
#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include "common.h"
#include "model.h"

class Controller {
public:

    Controller(Model::state s, Model::state targets);
    ~Controller() {
    }
    typedef std::vector<int> control;

    double calAssignment(Model::state s, Model::state targets);
    void calControl2d(Model::state s, Model::state targets);
    void translate_2d(double phi,Model::state s);
    double getDeviation(){return deviation;}

private:

    int  numTargets;
    double deviation;
    int dimP, numP;
    double radius;
    void calEudDist(Model::state s);
    std::vector<int> assignment;
    std::vector<double> availControl;
    Model::state targets_, s_;
    void calInlier(Model::state s);
};
