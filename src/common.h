#pragma once
#include <string>

struct Parameter{
	int N, dim, trajOutputInterval;
	double radius, dt, diffu_t, diffu_r, Bpp, Os_pressure, L_dep, cutoff,kappa;
	int captureStep, nCycles;
	std::string iniConfig, filetag, targetConfig;
    int seed;
    int captureOnlyFlag, captureAndTransportFlag; 
    int collective_MoveStep, collective_RestoreStep, CollectiveMoveCycle;
    int cargoCaptureStep;
    int assignViaEud;
    int dynamicTargetFlag;
    double controlTimeInterval, targetDiffuseRatio,  targetVelocityRatio, assignmentTimeInterval;
    double maxVelocity;
    int cargoInteractingFlag;
    double targetCenter[3];
    double totalCost;
        
    // transporter selection criterion
    int transporter_nb_thresh;
    double transporter_dist_thresh, transport_angle_thresh; // angle are in unit of degree
        
        
};

