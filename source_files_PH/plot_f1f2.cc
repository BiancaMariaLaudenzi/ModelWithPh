#include <fstream>
#include <iostream>
#include <tuple>
#include <iomanip>
#include "transport.h"
gasTransport gas; 

int main() {
    std::ofstream out("f1f2_vs_pO2pCO2_newton1.csv");
    out << "pO2,pCO2\n";

    double pCO2_min  = 30 * 0.133322;;  
    double pCO2_max  = 100 * 0.133322;;    
    double pO2_min   = 30 * 0.133322;;  
    double pO2_max   = 140 * 0.133322;;    
    double H_old     = 0.0000000339;  // mol/l
    double cCO2_target =  0.535371 / 22.7; // 0.0237; // 0.53799; //
    double cO2_target = 0.212045 / 22.7; // 0.009258; ; //0.2101566
    double pCO2_target = 38.9974 + 1; // 40; //* 0.133322;  // mmHg
    double pO2_target =  104.082 + 4; // 60; //* 0.133322;  // mmHg
    int    n_samples = 200;
    double cO2, cCO2;
    double pO2 =  pO2_target * 0.133322; 
	double pCO2 =  pCO2_target * 0.133322;
    cout << pO2 << " " << pCO2 << endl;

    const double tol = 1e-4;
    const double hCO2 = 1e-4 * pCO2;
	const double hO2 = 1e-4 * pO2;
    double cCO2_old = cCO2_target + 0.1 * cCO2_target;
	// convertire pressioni da mmHg a kPa 

    out << pO2 << "," << pCO2 << "\n";

    for (int iter = 0; iter < 100; ++iter) {
        
        double cCO2_pluspO2, cCO2_minuspO2, cO2_pluspCO2, cO2_minuspCO2, cO2_minuspO2, cO2_pluspO2, cCO2_minuspCO2, cCO2_pluspCO2;

    	tie(cO2, cCO2) = gas.dissociation(pCO2, pO2, cCO2_old, H_old);
		tie(cO2_pluspO2, cCO2_pluspO2) = gas.dissociation(pCO2, pO2 + hO2, cCO2_old, H_old);
		tie(cO2_minuspO2, cCO2_minuspO2) = gas.dissociation(pCO2, pO2 - hO2, cCO2_old, H_old);
		tie(cO2_pluspCO2, cCO2_pluspCO2) = gas.dissociation(pCO2 + hCO2, pO2, cCO2_old, H_old);
		tie(cO2_minuspCO2, cCO2_minuspCO2) = gas.dissociation(pCO2 - hCO2, pO2, cCO2_old, H_old);

		double f1 = (cO2 - cO2_target); //cO2_target;
        double f2 = (cCO2 - cCO2_target); //cCO2_target;

        // Derivate parziali numeriche
        double df1_pO2   = (cO2_pluspO2 - cO2_minuspO2) / (2*hO2);
        double df1_pCO2  = (cO2_pluspCO2 - cO2_minuspCO2) / (2*hCO2);
        double df2_pO2  = (cCO2_pluspO2 - cCO2_minuspO2) / (2*hO2);
        double df2_pCO2 = (cCO2_pluspCO2 - cCO2_minuspCO2) / (2*hCO2);
		// cout << "f1 = " << f1 << ", f2 = " << f2 << endl;
		// cout << "df1_pO2 = " << df1_pO2 << ", df1_pCO2 = " << df1_pCO2 << endl;
		// cout << "df2_pO2 = " << df2_pO2 << ", df2_pCO2 = " << df2_pCO2 << endl;

        // Matrice Jacobiana
        double det = df2_pCO2 * df1_pO2 - df2_pO2 * df1_pCO2;
        if (fabs(det) < 1e-5){
			// cout<<"Determinant is too small, division by zero avoided."<<endl;
		}; 
		// cout << "det = " << det << endl;

        // Inversa della Jacobiana * [f1, f2]
        double dpO2  = (f1 * df2_pCO2 / det - f2 * df1_pCO2 / det);
        double dpCO2 = (-f1 * df2_pO2 / det + f2 * df1_pO2 / det);

		// //cout << "dpO2 = " << dpO2 << ", dpCO2 = " << dpCO2 << endl;

        double alpha = 1.0;
		// for (int j = 0; j < 100; ++j) {
		// 	double new_pO2 = pO2 - alpha * dpO2;
		// 	double new_pCO2 = pCO2 - alpha * dpCO2;
		// 	double new_cO2, new_cCO2;

		// 	tie(new_cO2, new_cCO2) = gas.dissociation(new_pCO2, new_pO2, cCO2_old, H_old);
		// 	double f1_new = (new_cO2 - cO2_target); ///cO2_target;
		// 	double f2_new = (new_cCO2 - cCO2_target); ///cCO2_target;

		// 	double norm_f_new = std::max(std::abs(f1_new), std::abs(f2_new));
		// 	double norm_f_old = std::max(std::abs(f1), std::abs(f2));

		// 	if (norm_f_new < norm_f_old) break; // range pressioni
			
		// 	alpha *= 0.8;
		// }

		pCO2 -= alpha * dpCO2;
		pO2  -= alpha * dpO2;

        out << pO2 << "," << pCO2 << "\n";

        double rel_err = std::max(abs(f1) / cO2_target, abs(f2) / cCO2_target);
        if (rel_err < tol)
            break;
    }

	// convertire pressioni da kPa a mmHg
	double pO2_out =  pO2 / 0.133322; 
	double pCO2_out =  pCO2 / 0.133322;
	// cout << "pO2_out = " << pO2_out << ", pCO2_out = " << pCO2_out << endl;








    // tie(cO2, cCO2) = gas.dissociation(pCO2_target, pO2_target, cCO2_target+1e-2, H_old);
    // cout << "cO2 = " << cO2 << ", cCO2 = " << cCO2 << endl;

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

    // for (int i = 0; i <= n_samples; ++i) {
    //     for (int j = 0; j <= n_samples; ++j) {
    //         double pO2 = pO2_min + (pO2_max - pO2_min) * i / double(n_samples);
    //         double pCO2 = pCO2_min + (pCO2_max - pCO2_min) * j / double(n_samples);
    //         tie(cO2, cCO2) = gas.dissociation(pCO2, pO2, cCO2_target+1e-2, H_old);
    //         // cout << "cO2 = " << cO2 << ", cCO2 = " << cCO2 << endl;

    //         double f1 = cO2 - cO2_target;
    //         double f2 = cCO2 - cCO2_target;
    //         out << pO2 << "," << pCO2 << "," << f1 << "," << f2 << "\n";
    //     }
    // }
    // out.close();

    return 0;
}
