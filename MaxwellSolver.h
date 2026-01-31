#ifndef MaxwellSolver_H
#define MaxwellSolver_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <omp.h>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable> 
#include <atomic>



struct ParametersMaxwell{
    //Maxwell
    int nX, nY, nZ;
    double dt;
    double Delta;
    double epsilon;
    double mi;
    double q;
    double madderCoef;
    std::string B_OUTPUT;
    std::string E_OUTPUT;
    std::string pathTheta;
    std::string pathOmega;

    //pendulum
    double pTheta0;
    double pOmega0;
    double pL;
    double pMass;
    double pG;
    double tMax;

    //SIMULATIONPARAMETERS
    int NUM_OF_THREADS;
    };
struct field{
    double x;
    double y;
    double z;
};

struct SaveTask {
    int timeStep;
    std::vector<field> dataE;
    std::vector<field> dataB;
    std::string pathE;
    std::string pathB;
};



    class Solver{

        private:
            std::vector<field> m_E; // l= i + nX* j + nX*nY* k
            std::vector<field> m_E_old; // l= i + nX* j + nX*nY* k
            std::vector<field> m_B; // l= i + nX* j + nX*nY* k
            std::vector<field> m_currentDensity; // l= i + nX* j + nX*nY* k
            std::vector<field> m_PendullumPosition; //vector n*n*n
            std::vector<field> m_PendullumVelocity; //vector n*n*n
            std::vector<double> m_rho;
            std::vector<double> m_error;

            
            ParametersMaxwell m_parameters;

            std::ofstream m_file_E;
            std::ofstream m_file_B;

                //saver
            std::queue<SaveTask> m_saveQueue;
            std::mutex m_queueMtx;
            std::condition_variable m_cv;
            std::atomic<bool> m_running{true};
            std::thread m_writerThread;
            



// ============================== converts (i,j,k) to l ============================
            int converTol(int i ,int j, int k) const; 
// =============================== SAVING DATA WITHOUT STOPING CALC ===========
            void backgroundWriter();
        
        
            public:
// ============================== CONSTRUCTOR ======================================
            Solver(ParametersMaxwell parameters, const std::vector<field> &r, const std::vector<field> &v);

// =============================== DESTRUCTOR ========================================
            ~Solver();
// ============================== BLURS CURRENT DENSITY AND POINT CHARGE ON NEIGHBOURS ===========
            void Blur(int timeStep);
// ============================== RELAXATION IN SINGLE TIME STEP =======================================
            void relaxation(int timeStep);
// ============================== SAVES E AND B ==========================================
            void frontgroundWriter(int timeStep, int saveFreq);

    };
#endif

//README

//napisac praser do czytania plikow inputowych

//napisac interpolacje do zmiany dt wzgledem wahadla

