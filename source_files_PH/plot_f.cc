#include <fstream>
#include <iostream>
#include <tuple>
#include "transport.h"
gasTransport gas; 

int main() {
    // std::ofstream out("f_vs_cCO2.csv");
    // out << "cCO2,f\n";

    // double pCO2      = 40 * 0.133322; // kPa
    // double pO2       = 46 * 0.133322;  // kPa
    // double H_old     = 0.0000000339;  // mol/l
    // double cCO2_min  = 0.01;  // mol/l
    // double cCO2_max  = 0.03;    // mol/l
    // int    n_samples = 200;

    // for (int k = 0; k <= n_samples; ++k) {
    //     double cCO2 = cCO2_min + (cCO2_max - cCO2_min) * k / double(n_samples);
    //     double pH    = get<0>(gas.calculate_pH(cCO2, H_old));
    //     cout << "pH = " << pH << endl;
    //     double SO2   = gas.calculate_SO2(pO2, pCO2, pH);
    //     double phi   = gas.calculate_cCO2(pCO2, pH, SO2);
    //     double f     = phi - cCO2;
    //     out << cCO2 << "," << f << "\n";
    // }
    // out.close();

    // return 0;



    const double tol = 1e-4 * 0.022;
    const int max_iter = 100;
    double pCO2      = 40 * 0.133322; // kPa
    double pO2       = 46 * 0.133322;  // kPa
    double H_old     = 0.0000000339;  // mol/l

    

    double cCO2 = 0.03; // 22.0e-3 mol/l = 0.4994 mLgas/mLblood;  // iniziale
    std::ofstream out("f_vs_cCO2_NEWTON3.csv");
    out << "cCO2,f\n";
    for (int i = 0; i < max_iter; ++i) {
        double pH = get<0>(gas.calculate_pH(cCO2, H_old));
        double SO2 = gas.calculate_SO2(pO2, pCO2, pH);
        double phi = gas.calculate_cCO2(pCO2, pH, SO2);
        double f = phi - cCO2;
        out << cCO2 << "," << f << "\n";
        double df = gas.df_dcCO2(cCO2, pCO2, SO2, H_old);
        
        // Newton-Raphson step
        double delta = f / df;
        cCO2 -= delta;

        if (fabs(delta) < tol) break;
    }

}