#include "model.h"
#include "common.h"
#include <algorithm>
#define M_PI 3.1415926
double const Model::T = 293.0;
double const Model::kb = 1.38e-23;
double const Model::vis = 1e-3;


extern Parameter parameter;

Model::Model(){
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);

   
    filetag = parameter.filetag;
    iniFile = parameter.iniConfig;
    numP = parameter.N;
    dimP = 2;
    radius = parameter.radius;
    dt_ = parameter.dt;
    diffusivity_t = parameter.diffu_t;// this corresponds the diffusivity of 1um particle
    diffusivity_r = parameter.diffu_r; // this correponds to rotation diffusity of 1um particle
    Bpp = parameter.Bpp * kb * T * 1e9; //2.29 is Bpp/a/kT
    Kappa = parameter.kappa; // here is kappa*radius
    Os_pressure = parameter.Os_pressure * kb * T * 1e9;
    Os_pressure_origin = parameter.Os_pressure * kb * T * 1e9;
    L_dep = parameter.L_dep; // 0.2 of radius size, i.e. 200 nm
    radius_nm = radius*1e9;
    combinedSize = (1+L_dep)*radius_nm;
    mobility = diffusivity_t/kb/T;
    trajOutputInterval = parameter.trajOutputInterval;
    fileCounter = 0;
    cutoff = parameter.cutoff;
    this->rand_generator.seed(parameter.seed);

    for(int i = 0; i < numP; i++){
        particles.push_back(particle_ptr(new Model::particle));
        targets.push_back(particle_ptr(new Model::particle));
    }
       
}

void Model::run() {
    if (this->timeCounter == 0 || ((this->timeCounter + 1) % trajOutputInterval == 0)) {       
        this->outputTrajectory(this->trajOs);
    }
    
    calForces();

    for (int i = 0; i < numP; i++) {

        double random_x = sqrt(2.0 * diffusivity_t / dt_) * (*rand_normal)(rand_generator);


        double random_y = sqrt(2.0 * diffusivity_t / dt_) * (*rand_normal)(rand_generator);


        double ux = mobility * particles[i]->F[0] +
            parameter.maxVelocity * particles[i]->u * cos(particles[i]->phi)
            + random_x;


        double uy = mobility * particles[i]->F[1] +
            parameter.maxVelocity * particles[i]->u * sin(particles[i]->phi)
            + random_y;


        particles[i]->r[0] += ux * dt_;
        particles[i]->r[1] += uy * dt_;

        particles[i]->Fx = mobility * particles[i]->F[0];
        particles[i]->Fy = mobility * particles[i]->F[1];

        particles[i]->Vx = ux;
        particles[i]->Vy = uy;

        particles[i]->phi += sqrt(2.0 * diffusivity_r * dt_) * (*rand_normal)(rand_generator);
    }
    if (parameter.dynamicTargetFlag) {

        double randDist[2], disp[2];
        randDist[0] = (*rand_normal)(rand_generator);
        randDist[1] = (*rand_normal)(rand_generator);

        disp[0] = parameter.targetDiffuseRatio * mobility * targetCenter.F[0] * dt_
            + (parameter.targetVelocityRatio * parameter.maxVelocity * dt_ +
                sqrt(2.0 * diffusivity_t * parameter.targetDiffuseRatio * dt_) * randDist[0]);
        disp[1] = parameter.targetDiffuseRatio * mobility * targetCenter.F[1] * dt_
            + sqrt(2.0 * diffusivity_t * parameter.targetDiffuseRatio * dt_) * randDist[1];

        // change targetCenter; disp is in the unit of m
        targetCenter.r[0] += disp[0];
        targetCenter.r[1] += disp[1];

        // shift all targets 
        for (int i = 0; i < numP; i++) {
            targets[i]->r[0] += disp[0] / radius;
            targets[i]->r[1] += disp[1] / radius;
        }
    }
    this->timeCounter++;
}

void Model::run(int steps){
        
    for (int i = 0; i < steps; i++){
	    run();
    }
}



// this force calculation includes double layer repulsion and depletion attraction 
void Model::calForcesHelper_DLAO(double ri[3], double rj[3], double F[3],int i,int j,double& pot) {
    double r[2], dist;

    dist = 0.0;
    for (int k = 0; k < dimP; k++) {
        F[k] = 0.0;
        r[k] = (rj[k] - ri[k]) / radius;
        dist += pow(r[k], 2.0);
    }
    dist = sqrt(dist);
    if (dist < 2.0) {
        std::cerr << "overlap " << i << "\t" << j << "\t"<< this->timeCounter << "dist: " << dist <<std::endl;
        dist = 2.06;
    }
    if (dist < cutoff) {
        double Fpp = -4.0/3.0*
        Os_pressure* M_PI *(-3.0/4.0* pow(combinedSize,2.0)+3.0*dist*dist/16.0*radius_nm*radius_nm);
        Fpp += -Bpp * Kappa * exp(-Kappa*(dist-2.0));
        for (int k = 0; k < dimP; k++) {
            F[k] = Fpp * r[k] / dist;

        }
        //    Bpp = parameter.Bpp * kb * T * 1e9; //2.29 is Bpp/a/kT, a in unit of nm
        pot = parameter.Bpp * radius_nm * exp(-Kappa*(dist-2.0));
        double norm_comb = 1 +  L_dep;
        double ao = -4.0/3.0*parameter.Os_pressure*pow(radius_nm,3) *M_PI*(pow(norm_comb,3) -3.0/4.0*pow(norm_comb,2) *dist+1.0/16.0* pow(dist,3));
        pot += ao;
    }
}

// this force calculation only includes double layer repulsion 
void Model::calForcesHelper_DL(double ri[3], double rj[3], double F[3],int i, int j) {
    double r[2], dist;

    dist = 0.0;
    for (int k = 0; k < dimP; k++) {
        F[k] = 0.0;
        r[k] = (rj[k] - ri[k]) / radius;
        dist += pow(r[k], 2.0);
    }
    dist = sqrt(dist);
    if (dist < 2.0) {
        std::cerr << "overlap " << i << "\t with " << j << "\t"<< this->timeCounter << "dist: " << dist <<std::endl;
        dist = 2.06;
    }
    if (dist < cutoff) {
        double Fpp = -Bpp * Kappa * exp(-Kappa*(dist-1.9));
        
        for (int k = 0; k < dimP; k++) {
            F[k] = Fpp * r[k] / dist;
        }
    }
}

void Model::calForces() {
    double r[3], dist, F[3];
    for (int i = 0; i < numP; i++) {
        for (int k = 0; k < dimP; k++) {
            particles[i]->F[k] = 0.0;
        }
    }
    

    for (int i = 0; i < numP - 1; i++) {
        for (int j = i + 1; j < numP; j++) {
            double pot;
            calForcesHelper_DLAO(particles[i]->r, particles[j]->r, F,i, j,pot);
               
            for (int k = 0; k < dimP; k++) {
                particles[i]->F[k] += F[k];
                particles[j]->F[k] += -F[k];
            }
        }
    }
        
    if (parameter.cargoInteractingFlag){
        // Note here we treat the target center as the cargo
        for (int k = 0; k < dimP; k++) {
            targetCenter.F[k] = 0.0;
        }
            
        for (int j = 0; j < numP; j++) {
            calForcesHelper_DL(targetCenter.r, particles[j]->r, F,-1, j);
            for (int k = 0; k < dimP; k++) {
                targetCenter.F[k] += F[k];
                particles[j]->F[k] += -F[k];
            }
        }
        
    }

}
    


void Model::createInitialState(){

    this->readxyz(iniFile);
    this->readTarget(parameter.targetConfig);

    // target center is just the cargo position
    this->targetCenter.r[0] = parameter.targetCenter[0] * radius;
    this->targetCenter.r[1] = parameter.targetCenter[1] * radius;
    this->targetCenter.r[2] = parameter.targetCenter[2] * radius;
    
    std::stringstream ss;
    std::cout << "model initialize at round " << fileCounter << std::endl;
    ss << this->fileCounter++;
    if (trajOs.is_open()) trajOs.close();
    if (osCargo.is_open()) osCargo.close();
    this->trajOs.open(filetag + "xyz_" + ss.str() + ".txt");
    this->osTarget.open(filetag +"target"+ss.str() + ".txt");
    this->osCargo.open(filetag +"cargo"+ss.str() + ".txt");
    this->timeCounter = 0;
}

void Model::outputTrajectory(std::ostream& os) 
{

    for (int i = 0; i < numP; i++) {
        os << i << "\t";
        this->osTarget << i << "\t";
        for (int j = 0; j < 3; j++){
            os << particles[i]->r[j]/radius << "\t";
        }
        // first four columns id, x, y z
        for (int j = 0; j < 3; j++){
            osTarget << targets[i]->r[j] << "\t";
        }
        osTarget << std::endl;
        
        os << particles[i]->phi<< "\t";
        os << particles[i]->theta<< "\t";
        os << particles[i]->cost<< "\t";
        os << particles[i]->u<< "\t";
        os << particles[i]->targetIdx<< "\t";
        os << particles[i]->EudDistToTarget<< "\t";
        os << targets[i]->r[0]<< "\t";
        os << targets[i]->r[1]<< "\t";
        os << targets[i]->r[2]<<"\t";
        os << this->timeCounter*this->dt_ << "\t";
        os << particles[i]->transporterFlag<<"\t";

        os << std::endl;
    }
    for (int j = 0; j < 3; j++){
        this->osCargo << targetCenter.r[j]/radius << "\t";
    }
    this->osCargo << this->timeCounter*this->dt_ << "\t";
    this->osCargo << std::endl;
    
}



void Model::readxyz(const std::string filename) {
    std::ifstream is;
    is.open(filename.c_str());
    std::string line;
    double dum;
    for (int i = 0; i < numP; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> particles[i]->r[0];
        linestream >> particles[i]->r[1];
        linestream >> particles[i]->r[2];
        linestream >> particles[i]->phi;
        linestream >> particles[i]->theta;
    }
    for (int i = 0; i < numP; i++) {
        particles[i]->r[0] *=radius;
        particles[i]->r[1] *=radius;
        particles[i]->r[2] *=radius;
      
    }
    
    is.close();
}

void Model::readTarget(std::string filename) {
    std::ifstream is;
    is.open(filename);
    std::string line;
    std::stringstream linestream;
    double dum;
    for (int i = 0; i < numP; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> targets[i]->r[0];
        linestream >> targets[i]->r[1];
        linestream >> targets[i]->r[2];
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        // now do the target shift respect to the target center
        targets[i]->r[0] += parameter.targetCenter[0];
        targets[i]->r[1] += parameter.targetCenter[1];
        targets[i]->r[2] += parameter.targetCenter[2];

    }

}