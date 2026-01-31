#include "Interpol.h"

LagrangeInterpolation::LagrangeInterpolation(int sampleStart, int sampleEnd, std::vector<double> x,std::vector<double>function): 
                        m_sampleStart(sampleStart), m_sampleEnd(sampleEnd)
                        {
                            m_f.assign(function.begin() + sampleStart, function.begin() + sampleEnd + 1);
                            m_x.assign(x.begin()+sampleStart, x.begin()+sampleEnd+1);
                        }
                    

LagrangeInterpolation::~LagrangeInterpolation(){}

double LagrangeInterpolation::L_function(double x)
{
    int iterationNumb=m_sampleStart-m_sampleEnd+1;
    double Lf=0;
    double w=1.0;

    for (int i=0; i<iterationNumb; i++)
    {
        for (int j=0; j<iterationNumb; j++)
        {
            if (i==j) 
                continue;

            w*=(x-m_x[j])/(m_x[i]-m_x[j]);
        }
        w*=m_f[i];
        Lf+=w;
        w=1.0;
    }
    return Lf;
}