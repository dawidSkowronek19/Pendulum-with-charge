#include "pendulumEng.h"


PendulumEngine::PendulumEngine(ParametersMaxwell p)
    : m_parameters(p) {

    m_file_theta.open(p.pathTheta);
    m_file_omega.open(p.pathOmega);
    //m_file_Energy.open("/Energy.dat");

    m_theta = p.pTheta0;
    m_omega = p.pOmega0;
    m_file_theta<<std::fixed<<std::setprecision(12);
    m_file_omega<<std::fixed<<std::setprecision(12);
}
PendulumEngine::~PendulumEngine(){
    m_file_theta.close();
    m_file_omega.close();
    m_file_Energy.close();
}

//======================== Saving t, r(t), v(t) ============================/
void PendulumEngine::saveState(double theta, double omega, double t)
{    
    m_file_theta<< t << " " << theta << "\n";
    m_file_omega<< t << " " << omega << "\n";
    
}
//=============================== Saving K(t), V(t) ==============================/
/*void PendulumEngine::saveEnergy(double theta, double omega, double t)
{
    double K = m_parameters.m * m_parameters.L * m_parameters.L * omega * omega / 2.0;
    double V = -m_parameters.m * m_parameters.L *m_parameters.g * cos(theta);
    m_file_Energy << t << " " << K << " " << V <<std::endl;
}
*/
//======================= Velocity Verlet =================================/


void PendulumEngine::velocityVerlet()
{
    int iterrationNumber=m_parameters.tMax/m_parameters.dt;
    double force, forceOld;
    double t=0.0;
    forceOld=- m_parameters.pG * sin(m_theta)/m_parameters.pL;
    
    for (int step =0; step <iterrationNumber; step ++)
    {
        m_theta+=m_parameters.dt*m_omega + m_parameters.dt * m_parameters.dt * forceOld/2.0; 
        force=-m_parameters.pG * sin(m_theta)/m_parameters.pL;
        m_omega+= m_parameters.dt * (force + forceOld)/2.0;
        forceOld=force;
        t+=m_parameters.dt;
        saveState(m_theta,m_omega,t);

    }
}