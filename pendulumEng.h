#ifndef PENDULUM_H
#define PENDULUM_H

#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <cmath>
#include <iomanip>
#include "MaxwellSolver.h"

class PendulumEngine{
    public:
        PendulumEngine(ParametersMaxwell p);
        ~PendulumEngine();

        void velocityVerlet();
    private:
        double m_theta;
        double m_omega;
        ParametersMaxwell m_parameters;

        std::ofstream m_file_theta;
        std::ofstream m_file_omega;
        std::ofstream m_file_Energy;
//======================== Saving t, r(t), v(t) ============================/
        void saveState(double theta, double omega, double t);

//======================== Save Energy ===================================/
        //void saveEnergy(double theta, double omega, double t);

};

#endif