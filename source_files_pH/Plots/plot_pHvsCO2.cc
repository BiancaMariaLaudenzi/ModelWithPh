#include <fstream>
#include <iostream>
#include <tuple>
#include <iomanip>
#include "transport.h"
gasTransport gas; 

int main() {
    std::ofstream out("pCO2vscCO2_pO25.3.csv");
    out << "pCO2,cCO2\n";

    double H_old     = 0.0000000339;  // mol/l
    double cCO2_min =  0.005; // 0.0237; 
    double cCO2_max = 0.05; // 0.009258; ;
    double pCO2_min = 1.0 * 0.133322; // kPa
    double pCO2_max =  100 * 0.133322; // kPa
    double cCO2_old = 0.0237; // mol/l
    int    n_samples = 200;
    double cO2, cCO2;

    double pH, H;

    for (int i = 0; i <= n_samples; ++i) {
        double pO2 = 5.3; // kPa
        double pCO2 = pCO2_min + (pCO2_max - pCO2_min) * i / double(n_samples);
        tie(cO2, cCO2) = gas.dissociation(pCO2, pO2, cCO2_old, H_old);

        out << pCO2 << "," << cCO2 << "\n";
    }
    out.close();

    // double pppo2 = 104.082, pppco2 = 38.9974;
    // double pepo2 = 67.5991, pepco2 = 37.6065;
    // double pspo2 = 30.6311, pspco2 = 45.6982;
    // double pmpo2 = 38.0859, pmpco2 = 40.7813;
    // double phpo2 = 40.5419, phpco2 = 49.6231;
    // double pbpo2 = 32.7514, pbpco2 = 43.9942;
    // double pevo2 = 67.5991, pevco2 = 37.6065;
    // double psvo2 = 30.6311, psvco2 = 45.6982;
    // double pmvo2 = 38.0859, pmvco2 = 40.7813;
    // double phvo2 = 40.5419, phvco2 = 49.6231;
    // double pbvo2 = 32.7514, pbvco2 = 43.9942;
    // double pvo2  = 39.0809, pvco2 = 42.9239;

    // double cppo2 = 0.212045, cppco2 = 0.535371;
    // double cepo2 = 0.210520, cepco2 = 0.529065;
    // double cspo2 = 0.208371, cspco2 = 0.563870;
    // double cmpo2 = 0.208917, cmpco2 = 0.543290;
    // double chpo2 = 0.209078, chpco2 = 0.579201;
    // double cbpo2 = 0.208537, cbpco2 = 0.556919;
    // double cevo2 = 0.210520, cevco2 = 0.529065;
    // double csvo2 = 0.208371, csvco2 = 0.563870;
    // double cmvo2 = 0.208917, cmvco2 = 0.543290;
    // double chvo2 = 0.209078, chvco2 = 0.579201;
    // double cbvo2 = 0.208537, cbvco2 = 0.556919;
    // double cvo2  = 0.208983, cvco2  = 0.552452;

    // auto calc_print = [&](std::string nome, double cO2, double cCO2, double pO2, double pCO2) {
    //     std::tie(pCO2, pO2) = gas.invertedDissociation(cO2, cCO2, H_old, pO2 + 0.1 * pO2, pCO2 + 0.1 * pCO2);
    //     std::cout << nome << ": pO2 = " << pO2 << " mmHg, pCO2 = " << pCO2 << " mmHg" << std::endl;
    // };

    // // Chiamate per tutti i punti
    // calc_print("cpp", cppo2, cppco2, pppo2, pppco2);
    // calc_print("cep", cepo2, cepco2, pepo2, pepco2);
    // calc_print("csp", cspo2, cspco2, pspo2, pspco2);
    // calc_print("cmp", cmpo2, cmpco2, pmpo2, pmpco2);
    // calc_print("chp", chpo2, chpco2);
    // calc_print("cbp", cbpo2, cbpco2);
    // calc_print("cev", cevo2, cevco2);
    // calc_print("csv", csvo2, csvco2);
    // calc_print("cmv", cmvo2, cmvco2);
    // calc_print("chv", chvo2, chvco2);
    // calc_print("cbv", cbvo2, cbvco2);
    // calc_print("cv",  cvo2,  cvco2);

    return 0;
}
