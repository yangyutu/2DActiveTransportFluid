#pragma once
#include<vector>
#include<memory>
#include<random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>


class Model {
public:    
    struct particle {
        double r[3]{0.0, 0.0, 0.0}, F[3]{ 0.0, 0.0, 0.0 }; // positions and forces
        double phi{ 0.0 }; // angle phi
        double theta{ 0.0 }; // angle theta
        double u {0.0}; // normalized propulsion speed, within range [0, 1] 
        int targetIdx{-1}; // target to track
        double targetPos[3]{ 0.0, 0.0, 0.0 }; // target position to track
        double EudDistToTarget{ 0.0 }; // Euclidean distance to the target to be tracked
        double cost;
        int nbcount; // number of neighbors
        int inlier; // whether it is an inlier

        int transporterFlag;      // whether it is a transporter  


        double Fx{ 0.0 }, Fy{ 0.0 }, Vx{ 0.0 }, Vy{0.0}; // Fx, Fy are the drifting velocity
                
        
        particle(double x = 0, double y = 0, double z = 0){
            r[0]=x; r[1]=y; r[2]=z;
        }
    };
    typedef std::shared_ptr<particle> particle_ptr;
    typedef std::vector<particle_ptr> state;

   
    Model();
    ~Model() {
        trajOs.close();
        osTarget.close();
    }
    void run();
    void run(int steps);
    void createInitialState();
    state getCurrState(){return particles;}
    state getTargets(){return targets;}
    double dt(){return dt_;}
    int np(){return numP;}
    
private:
    void calForces();
    void calForcesHelper_DLAO(double ri[3], double rj[3], double F[3],int i, int j, double& pot);
    void calForcesHelper_DL(double ri[3], double rj[3], double F[3],int i, int j);
    

    // physical parameters
    int dimP{2};
    static const double kb, T, vis;
    int numP;
    double radius, radius_nm;
    double LJ, rm;
    double Bpp; //2.29 is Bpp/a/kT
    double Kappa; // here is kappa*radius
    double Os_pressure, Os_pressure_origin;
    double L_dep; // 0.2 of radius size, i.e. 200 nm
    double combinedSize;
    double dt_, cutoff, mobility, diffusivity_r, diffusivity_t;

    // particle, target, and target center(cargo) positions
    state particles, targets,initialDistToCenter;
    particle targetCenter; // target center is the cargo position

    // random number
    std::default_random_engine rand_generator;
    std::shared_ptr<std::normal_distribution<double>> rand_normal;

    // data IO
    std::string iniFile;
    long long timeCounter, fileCounter;
    std::ofstream trajOs, osTarget, osCargo;
    std::string filetag;
    int trajOutputInterval;
    void outputTrajectory(std::ostream& os);
    void readxyz(const std::string filename);
    void readTarget(std::string filename);
    
};



