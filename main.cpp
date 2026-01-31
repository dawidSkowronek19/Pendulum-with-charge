#include "pendulumEng.h"
#include "MaxwellSolver.h"
#include "json.hpp"

using json = nlohmann::json;

class ConfigParser {
public:

    static ParametersMaxwell loadConfig(const std::string& filename) {
        std::ifstream file(filename);

        json j;
        file >> j;

        ParametersMaxwell p;
        p.nX = j["simulation"]["nX"];
        p.nY = j["simulation"]["nY"];
        p.nZ = j["simulation"]["nZ"];
        p.Delta = j["simulation"]["Delta"];
        p.dt = j["simulation"]["timeStep"];
        p.epsilon = j["simulation"]["epsilon"];
        p.mi = j ["simulation"]["mi"];
        p.q = j["simulation"]["q"];
        p.madderCoef = j["simulation"]["MadderCoef"];


        p.pTheta0 = j["pendulum"]["theta0"];
        p.pOmega0 = j["pendulum"]["omega0"];
        p.pL = j["pendulum"]["L"];
        p.pMass = j["pendulum"]["mass"];
        p.pG = j["pendulum"]["g"];
        p.tMax = j["pendulum"]["tMax"];


        p.E_OUTPUT = j["pathSimulation"]["E_OUTPUT"];
        p.B_OUTPUT = j["pathSimulation"]["B_OUTPUT"];
        p.pathTheta = j["pathPendulum"]["thetaFile"]; 
        p.pathOmega = j["pathPendulum"]["omegaFile"];

        p.NUM_OF_THREADS =j["simulationParameters"]["NUM_OF_THREADS"];

        return p;
    }

};






std::pair<std::vector<field>, std::vector<field>> loadFromSeparateFiles(ParametersMaxwell p) {
    std::vector<field> r, v;
    std::ifstream fileTheta(p.pathTheta);
    std::ifstream fileOmega(p.pathOmega);

    if (!fileTheta.is_open() || !fileOmega.is_open()) {
        throw std::runtime_error("FILES NOT FOUND\n");
    }

    double t_dummy, theta, omega;
    double offX = p.nX*p.Delta/2.0, offY = p.nY*p.Delta/2.0, offZ = p.nZ*p.Delta/2.0; 

    while (fileTheta >> t_dummy >> theta && fileOmega >> t_dummy >> omega) {
        double px = offX + p.pL * std::sin(theta);
        double py = offY - p.pL * std::cos(theta); 
        double pz = offZ;
        r.push_back({px, py, pz});

        double vx = p.pL * omega * std::cos(theta);
        double vy = p.pL * omega * std::sin(theta);
        double vz = 0.0;
        v.push_back({vx, vy, vz});
    }

    return {r, v};
}

int main()
{
    
    ParametersMaxwell TEST = ConfigParser::loadConfig("config.json");
    omp_set_num_threads(TEST.NUM_OF_THREADS);

    double c = 1.0/sqrt(TEST.epsilon*TEST.mi);

    if (TEST.dt>=0.5*TEST.Delta/c) //Courant Limit
    {
        std::cout<<"Courant conditions is not satisified\n";
        return 1;
    }
    
    {
        PendulumEngine pendulumTest(TEST);
        pendulumTest.velocityVerlet();
    }

    auto [r_vec, v_vec] = loadFromSeparateFiles(TEST);

   
    int timeIterationMax=TEST.tMax/TEST.dt;

    Solver fieldTEST(TEST, r_vec, v_vec);

    for (int timeStep=0; timeStep<timeIterationMax; timeStep++)
    {
        fieldTEST.Blur(timeStep);
        fieldTEST.relaxation(timeStep);
        fieldTEST.frontgroundWriter(timeStep, 50);
    }


    return 0;
}