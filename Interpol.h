#ifndef Interpol_h
#define Interpol_h

#include <vector>
#include <iostream>

class LagrangeInterpolation{

    public:
        LagrangeInterpolation(int sampleStart, int sampleEnd, std::vector<double> m_x, std::vector<double>);
        ~LagrangeInterpolation();
        double L_function(double x);

    private:
        int m_sampleStart; //including m_sampleStart
        int m_sampleEnd; //including m_sampleEnd
        std::vector<double> m_f;
        std::vector<double> m_x;

};

#endif