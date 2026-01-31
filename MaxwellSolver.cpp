#include "MaxwellSolver.h"

Solver::Solver(ParametersMaxwell parameters, const std::vector<field>&r, const std::vector<field> &v): 
m_parameters(parameters), m_PendullumPosition(r), m_PendullumVelocity(v) 
{
    int volume = m_parameters.nX * m_parameters.nY *m_parameters.nZ;
    m_E.assign(volume, {0.0, 0.0, 0.0});
    m_E_old.assign(volume, {0.0, 0.0, 0.0});
    m_B.assign(volume, {0.0, 0.0, 0.0});
    m_currentDensity.assign(volume, {0.0, 0.0, 0.0});
    m_rho.assign(volume, 0.0);
    m_error.assign(volume, 0.0);

    m_writerThread = std::thread(&Solver::backgroundWriter, this);



}

Solver::~Solver()
{
    m_running = false;
    m_cv.notify_one();
    if (m_writerThread.joinable()) 
    {
        m_writerThread.join();
    }
}

int Solver::converTol(int i, int j, int k) const {
    return i + m_parameters.nX*j + m_parameters.nX*m_parameters.nY*k;
}

void Solver::Blur(int timeStep)
{
    std::fill(m_rho.begin(), m_rho.end(), 0.0);
    std::fill(m_currentDensity.begin(), m_currentDensity.end(), field{0.0, 0.0, 0.0});
    
    double Delta=m_parameters.Delta;
    double q = m_parameters.q;
    int nX = m_parameters.nX;
    int nY = m_parameters.nY;
    int nZ = m_parameters.nZ;
    
    double i_pendullum = m_PendullumPosition[timeStep].x/Delta;
    double j_pendullum = m_PendullumPosition[timeStep].y/Delta;
    double k_pendullum = m_PendullumPosition[timeStep].z/Delta;

    double exponent;
    double s =3.5; // softening
    double invs2=1.0/(s*s);
    double normalization = q / (pow(Delta, 3) * pow(s * sqrt(M_PI), 3));//normalization of cloud


    std::vector<double> e_k, e_j, e_i;
    
    e_k.resize(nZ);
    e_j.resize(nY);
    e_i.resize(nX);

    

    const int k_start=(int)(k_pendullum-3.0*s); //for pragma
    const int k_end=(int)(k_pendullum+3.0*s); //for pragma

    #pragma omp parallel for
    for (int k=k_start; k<k_end; k++)
    {
        if (k < 0 || k >= m_parameters.nZ) continue;
        double tmp_k=-(k-k_pendullum)*(k-k_pendullum)*invs2;
        e_k[k]=exp(tmp_k);
        for (int j=(int)(j_pendullum-3.0*s); j<(int)(j_pendullum+3.0*s); j++)
        {
            if (j < 0 || j >= m_parameters.nY) continue;
            double tmp_j=-(j-j_pendullum)*(j-j_pendullum)*invs2;
            e_j[j]=exp(tmp_j);
            for (int i=(int)(i_pendullum-3.0*s); i<(int)(i_pendullum+3.0*s); i++)
            {
                if (i < 0 || i >= m_parameters.nX) continue;
                double tmp_i=-(i-i_pendullum)*(i-i_pendullum)*invs2;
                e_i[i]=exp(tmp_i);
                
                int l=converTol(i,j,k);
                exponent = e_k[k]*e_j[j]*e_i[i];
                m_currentDensity[l].x=m_PendullumVelocity[timeStep].x*exponent*normalization;
                m_currentDensity[l].y=m_PendullumVelocity[timeStep].y*exponent*normalization;
                m_currentDensity[l].z=0.0;

                m_rho[l]=exponent*normalization;
            }
        }
    }

}

// ======================== RELAXATION IN TIME STEP =====================
/*
--------------------------------------------------------------------------
It is assumed that time step from pendullum calculations is the same as 
the one set in parameters.

Possible to add interpolation later
-----------------------------------------------------------------------------

*/

void Solver::relaxation(int timeStep)
{
    double dt=m_parameters.dt;
    double epsilon=m_parameters.epsilon;
    double mi=m_parameters.mi;
    double q=m_parameters.q;
    double Delta=m_parameters.Delta;



    int nX = m_parameters.nX;
    int nY = m_parameters.nY;
    int nZ =m_parameters.nZ;

    double multiplicator=dt/Delta;


    /*
    ----------------------------------------------------------------------------------
    The grid in 1D is shown below

    E_0(wall) ----- B_0 ----- E_1 ----- B_1 ----- E_2 ------ B_2 ----- E_3(wall)
    
    to calc B we need to use forward difference, to calc E we need to use backward 
    difference in order to achieve symmetry
    -------------------------------------------------------------------------------
    */


    // we start from n-1/2 time step. now we calculate B in (n+1/2) time step from E in n time step
    #pragma omp parallel for collapse(3)
    for (int k=0; k<nZ-1; k++)
    {
        for (int j=0; j<nY-1; j++)
        {
            for (int i=0; i<nX-1; i++)
            {
                int l = converTol(i,j,k);
                m_B[l].x+=-multiplicator*(m_E[l+nX].z-m_E[l].z-(m_E[l+nX*nY].y-m_E[l].y));
                m_B[l].y+=-multiplicator*(m_E[l+nX*nY].x-m_E[l].x-(m_E[l+1].z-m_E[l].z));
                m_B[l].z+=-multiplicator*(m_E[l+1].y-m_E[l].y-(m_E[l+nX].x-m_E[l].x));

            }
        }
    }
    // now we need to calculate E in (n+1) time step from B in (n+1/2) time step and E in n time step
    
    
    #pragma omp parallel for
    for (int i = 0; i < (int)m_E.size(); i++) m_E_old[i] = m_E[i];
    
    multiplicator=dt/(mi*epsilon*Delta);
    #pragma omp parallel for collapse(3)
    for (int k=1; k<nZ-1; k++)
    {
        for (int j=1; j<nY-1; j++)
        {
            for (int i=1; i<nX-1; i++)
            {
                int l = converTol(i,j,k);

                m_E[l].x+=multiplicator*(m_B[l].z-m_B[l-nX].z-(m_B[l].y-m_B[l-nX*nY].y)-mi*Delta*m_currentDensity[l].x);
                m_E[l].y+=multiplicator*(m_B[l].x-m_B[l-nX*nY].x-(m_B[l].z-m_B[l-1].z)-mi*Delta*m_currentDensity[l].y);
                m_E[l].z+=multiplicator*(m_B[l].y-m_B[l-1].y-(m_B[l].x-m_B[l-nX].x)-mi*Delta*m_currentDensity[l].z);
            }
        }
    }

    //Madder correction
    #pragma omp parallel for collapse(3)
    for (int k = 1; k < nZ - 1; k++) 
    {
        for (int j = 1; j < nY - 1; j++) 
        {
            for (int i = 1; i < nX - 1; i++) 
            {
                int l = converTol(i,j,k);
                double divE = (m_E[l].x - m_E[l - 1].x) / Delta +
                              (m_E[l].y - m_E[l - nX].y) / Delta +
                              (m_E[l].z - m_E[l - nX * nY].z) / Delta;
                m_error[l] = divE - (m_rho[l] / epsilon);
            }
        }
    }

    double madderCoef= m_parameters.madderCoef*m_parameters.dt; 
    #pragma omp parallel for collapse(3)
    for (int k = 1; k < nZ - 1; k++) 
    {
        for (int j = 1; j < nY - 1; j++) 
        {
            for (int i = 1; i < nX - 1; i++) 
            {
                int l = i + nX * j + nX * nY * k;
                m_E[l].x += madderCoef * (m_error[l + 1] - m_error[l]);
                m_E[l].y += madderCoef * (m_error[l + nX] - m_error[l]);
                m_E[l].z += madderCoef * (m_error[l + nX * nY] - m_error[l]);
            }
        }
    }

    //================================= BOUNDARY ==================================
    double c=1.0/sqrt(mi*epsilon);
    #pragma omp parallel for collapse(2)
    for (int k=0; k<nZ; k++)
    {
        for (int j=0; j<nY; j++)
        {
            int i=0;
            int l=converTol(i,j,k);
            m_E[l].y=m_E_old[l+1].y+(c*dt-Delta)/(c*dt+Delta)*(m_E[l+1].y-m_E_old[l].y);
            m_E[l].z=m_E_old[l+1].z+(c*dt-Delta)/(c*dt+Delta)*(m_E[l+1].z-m_E_old[l].z);

            i=nX-1;
            l =converTol(i,j,k);
            m_E[l].y=m_E_old[l-1].y+(c*dt-Delta)/(c*dt+Delta)*(m_E[l-1].y-m_E_old[l].y);
            m_E[l].z=m_E_old.at(l-1).z+(c*dt-Delta)/(c*dt+Delta)*(m_E[l-1].z-m_E_old[l].z);
        }
    }
    #pragma omp parallel for collapse(2)
    for (int k=0; k<nZ; k++)
    {
        for (int i=0; i<nX; i++)
        {
            int j=0;
            int l = converTol(i,j,k);
            m_E[l].x=m_E_old[l+nX].x+(c*dt-Delta)/(c*dt+Delta)*(m_E[l+nX].x-m_E_old[l].x);
            m_E[l].z=m_E_old[l+nX].z+(c*dt-Delta)/(c*dt+Delta)*(m_E[l+nX].z-m_E_old[l].z);

            j= nY-1;
            l=converTol(i,j,k);
            m_E[l].x=m_E_old[l-nX].x+(c*dt-Delta)/(c*dt+Delta)*(m_E[l-nX].x-m_E_old[l].x);
            m_E[l].z=m_E_old[l-nX].z+(c*dt-Delta)/(c*dt+Delta)*(m_E[l-nX].z-m_E_old[l].z);
        }
    }
    #pragma omp parallel for collapse(2)
    for (int j=0; j<nY; j++)
    {
        for (int i=0; i<nX; i++)
        {
            int k=0;
            int l = converTol(i,j,k);
            m_E[l].x=m_E_old[l+nX*nY].x+(c*dt-Delta)/(c*dt+Delta)*(m_E[l+nX*nY].x-m_E_old[l].x);
            m_E[l].y=m_E_old[l+nX*nY].y+(c*dt-Delta)/(c*dt+Delta)*(m_E[l+nX*nY].y-m_E_old[l].y);

            k=nZ-1;
            l=converTol(i,j,k);
            m_E[l].x=m_E_old[l-nX*nY].x+(c*dt-Delta)/(c*dt+Delta)*(m_E[l-nX*nY].x-m_E_old[l].x);
            m_E[l].y=m_E_old[l-nX*nY].y+(c*dt-Delta)/(c*dt+Delta)*(m_E[l-nX*nY].y-m_E_old[l].y);
        }
    }
}

// =============================== SAVES E & B IN RAM ==================

void Solver::frontgroundWriter(int timeStep, int saveFreq)
{
    if (timeStep % saveFreq == 0) 
    {
        SaveTask task;
        task.timeStep = timeStep;
        task.pathE = m_parameters.E_OUTPUT;
        task.pathB = m_parameters.B_OUTPUT;
        
        task.dataE = m_E; 
        task.dataB = m_B;

        std::lock_guard<std::mutex> lock(m_queueMtx);
        m_saveQueue.push(std::move(task));
        m_cv.notify_one();
    }
}

// ============================ SAVES E & B FROM RAM QUEUE ON DISK ================
void Solver::backgroundWriter() 
{
    while (m_running || !m_saveQueue.empty()) {
        std::unique_lock<std::mutex> lock(m_queueMtx);
        
        m_cv.wait(lock, [this] { return !m_saveQueue.empty() || !m_running; });

        if (!m_saveQueue.empty()) {
            SaveTask task = std::move(m_saveQueue.front());
            m_saveQueue.pop();
            lock.unlock();
            std::string fileNameE = task.pathE + std::to_string(task.timeStep) + ".bin";
            std::ofstream fe(fileNameE, std::ios::binary);
            fe.write(reinterpret_cast<const char*>(task.dataE.data()), task.dataE.size() * sizeof(field));
            fe.close();

            std::string fileNameB = task.pathB + std::to_string(task.timeStep) + ".bin";
            std::ofstream fb(fileNameB, std::ios::binary);
            fb.write(reinterpret_cast<const char*>(task.dataB.data()), task.dataB.size() * sizeof(field));
            fb.close();
        }
    }
}
