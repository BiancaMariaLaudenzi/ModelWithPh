/*
 * Copyright 2025 Bianca Maria Laudenzi, Caterina Dalmaso, Lucas Omar Muller
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include "GetPot.h"
#include <string>
#include <tuple>

#include "transport.h"

using namespace std;
using std::scientific;

void gasTransport::init(string ifile, string ifileC, string outDir) {

	// read parameter file
	if (verbose)
		cout << "Reading GAS params ..." << endl;
	if (ifile != "NOFILE") {
        if (ifileC !="NOFILE"){
		    readGASparameters(ifile,ifileC);
        }
	}
	if (verbose)
		cout << "Done reading GAS params" << endl;
		
	// intialiaze output file
	nameGAS = "GAS";
    string statenameStateVars = outDir + nameGAS + "_stateVars.txt";
	sampleGASstateVars.open(statenameStateVars.c_str());
    sampleGASstateVars << "# [0]: time;" << " "
					<< "[1]: fdo2;" << " "
                    << "[2]: fdco2;" << " "
                    << "[3]: fAo2;" << " "
                    << "[4]: fAco2;" << " "
                    << "[5]: cppo2;" << " "
                    << "[6]: cppco2;" << " "
                    << "[7]: cepo2;" << " "
                    << "[8]: cepco2;" << " "
                    << "[9]: cspo2;" << " "
                    << "[10]: cspco2;" << " "
                    << "[11]: cmpo2;" << " "
                    << "[12]: cmpco2;" << " "
                    << "[13]: chpo2;" << " " 
                    << "[14]: chpco2;" << " " 
                    << "[15]: cbpo2;" << " " 
                    << "[16]: cbpco2;" << " " 
					<< "[17]: cevo2;" << " " 
					<< "[18]: cevco2;" << " " 
                    << "[19]: csvo2;" << " " 
                    << "[20]: csvco2;" << " " 
                    << "[21]: cmvo2;" << " " 
				    << "[22]: cmvco2;" << " " 
                    << "[23]: chvo2;" << " " 
                    << "[24]: chvco2;" << " " 
                    << "[25]: cbvo2;" << " " 
                    << "[26]: cbvco2;" << " " 
                    << "[27]: cvo2;" << " " 
                    << "[28]: cvco2;" << " " 
                    << "\n";


	string statenameAlveoliLungsAlg = outDir + nameGAS + "_AlveoliLungsAlg.txt";
	sampleGASalveoliLungsAlg.open(statenameAlveoliLungsAlg.c_str());
	sampleGASalveoliLungsAlg << "# [0]: time;" << " "
                    << "[1]: uvdot;" << " "
                    << "[2]: umvdot;" << " "
                    << "[3]: pAco2;" << " "
                    << "[4]: pppo2;" << " "
                    << "[5]: pppco2;" << " "
                    << "[6]: cao2;" << " "
                    << "[7]: caco2;" << " "
                    << "[8]: pao2;" << " "
                    << "[9]: paco2;" << " "
                    << "[10]: cao2diss;" << " "
                    << "[11]: mdoto2;" << " "
                    << "[12]: mdotco2;" << " "
					<< "[13]: pAo2;" << " "
					<< "[14]: pHsysArt;" << " "
					<< "[15]: pHpulCap;" << " "
					<< "\n";


	string statenameTissueVenousGasAlg = outDir + nameGAS + "_TissueVenousGasAlg(delayed).txt";
	sampleGAStissueVenousGasAlg.open(statenameTissueVenousGasAlg.c_str());
	sampleGAStissueVenousGasAlg << "# [0]: time;" << " "
                    << "[1]: cao2tilde;" << " "
                    << "[2]: caco2tilde;" << " "
                    << "[3]: cvo2tilde;" << " "
                    << "[4]: cvco2tilde;" << " "
					<< "\n";
					
	InitialCond(ifile);
	
	if (verbose){
		cout << "gasTransport::init :: DONE " << endl;
	}

}


void gasTransport::readGASparameters(string _file, string _fileC) {
	if (verbose)
		cout << " -- gasTransport::readGASparameters() reading from "<< _file << " and " << _fileC << endl;

	GetPot ifileC(_fileC.c_str());  // common parameters are assignes in commonParams.dat
	GetPot ifile(_file.c_str());
		
	// [Lungs Gas]
	LungsGasParams[0] = ifile("LungsGas/fio2",0.);
	LungsGasParams[1] = ifile("LungsGas/fico2",0.);
	LungsGasParams[2] = ifile("LungsGas/k",0.);
	LungsGasParams[3] = ifile("LungsGas/patm",0.);
	LungsGasParams[4] = ifile("LungsGas/pws", 0.);
	LungsGasParams[5] = ifile("LungsGas/csato2",0.);
	LungsGasParams[6] = ifile("LungsGas/csatco2",0.);
    LungsGasParams[7] = ifile("LungsGas/h1",0.);
    LungsGasParams[8] = ifile("LungsGas/h2",0.);
    LungsGasParams[9] = ifile("LungsGas/alpha1",0.);
    LungsGasParams[10] = ifile("LungsGas/alpha2",0.);
    LungsGasParams[11] = ifile("LungsGas/beta1",0.);
    LungsGasParams[12] = ifile("LungsGas/beta2",0.);
    LungsGasParams[13] = ifile("LungsGas/k1",0.);
    LungsGasParams[14] = ifile("LungsGas/k2",0.);
	LungsGasParams[15] = ifileC("shunt/sh", 0.);
	LungsGasParams[16] = ifile("LungsGas/hgb", 0.);
	LungsGasParams[17] = ifile("LungsGas/kO2", 0.);
	LungsGasParams[18] = ifile("LungsGas/kCO2", 0.);


	// [Tissues Venous Gas] 
	// N.B. consumption rates have to be converted in mL/seconds.
	TissuesVenousGasParam[0] = ifile("TissuesVenousGas/vtep",0.);
	TissuesVenousGasParam[1] = ifile("TissuesVenousGas/mo2ep",0.)/60.;
    TissuesVenousGasParam[2] = ifile("TissuesVenousGas/mco2ep",0.)/60.;
	TissuesVenousGasParam[3] = ifile("TissuesVenousGas/vtsp", 0.);
    TissuesVenousGasParam[4] = ifile("TissuesVenousGas/mo2sp",0.)/60.;
    TissuesVenousGasParam[5] = ifile("TissuesVenousGas/mco2sp",0.)/60.;
    TissuesVenousGasParam[6] = ifile("TissuesVenousGas/vtmp",0.);
    TissuesVenousGasParam[7] = ifile("TissuesVenousGas/mo2mp",0.)/60.;
    TissuesVenousGasParam[8] = ifile("TissuesVenousGas/mco2mp",0.)/60.;
    TissuesVenousGasParam[9] = ifile("TissuesVenousGas/vthp",0.);
    TissuesVenousGasParam[10] = ifile("TissuesVenousGas/mo2hp",0.)/60.;
    TissuesVenousGasParam[11] = ifile("TissuesVenousGas/mco2hp",0.)/60.;
    TissuesVenousGasParam[12] = ifile("TissuesVenousGas/vtbp",0.);
    TissuesVenousGasParam[13] = ifile("TissuesVenousGas/mo2bp",0.)/60.;
    TissuesVenousGasParam[14] = ifile("TissuesVenousGas/mco2bp", 0.)/60.;
	
	// [DelayedGasParam] (in closedloop.cc)
	DelayedGasParam[0] = ifile("DelayedGasParam/tault",0.);
	DelayedGasParam[1] = ifile("DelayedGasParam/tauvl",0.);

	if (verbose)
		printGASparameters();
}

void gasTransport::InitialCond(string _file) {
	if (verbose)
		cout << " -- gasTransport::saveInitialCond() reading from "
		<< _file << endl;

	GetPot ifile(_file.c_str());
	// [initialConditions]
	stateVars[0] = ifile("initial_conditions/fdo2", 0.);
	stateVars[1] = ifile("initial_conditions/fdco2", 0.);
	stateVars[2] = ifile("initial_conditions/fAo2", 0.);
	stateVars[3] = ifile("initial_conditions/fAco2", 0.);
	if (initial_cond = 0) {
	stateVars[4] = ifile("initial_conditions/cppo2", 0.);
	stateVars[5] = ifile("initial_conditions/cppco2", 0.);
	stateVars[6] = ifile("initial_conditions/cepo2", 0.);
	stateVars[7] = ifile("initial_conditions/cepco2", 0.);
	stateVars[8] = ifile("initial_conditions/cspo2", 0.);
	stateVars[9] = ifile("initial_conditions/cspco2", 0.);
	stateVars[10] = ifile("initial_conditions/cmpo2", 0.);
	stateVars[11] = ifile("initial_conditions/cmpco2", 0.);
	stateVars[12] = ifile("initial_conditions/chpo2", 0.);
	stateVars[13] = ifile("initial_conditions/chpco2", 0.);
	stateVars[14] = ifile("initial_conditions/cbpo2", 0.);
	stateVars[15] = ifile("initial_conditions/cbpco2", 0.);
	stateVars[16] = ifile("initial_conditions/cevo2", 0.);
	stateVars[17] = ifile("initial_conditions/cevco2", 0.);
	stateVars[18] = ifile("initial_conditions/csvo2", 0.);
	stateVars[19] = ifile("initial_conditions/csvco2", 0.);
	stateVars[20] = ifile("initial_conditions/cmvo2", 0.);
	stateVars[21] = ifile("initial_conditions/cmvco2", 0.);
	stateVars[22] = ifile("initial_conditions/chvo2", 0.);
	stateVars[23] = ifile("initial_conditions/chvco2", 0.);
	stateVars[24] = ifile("initial_conditions/cbvo2", 0.);
	stateVars[25] = ifile("initial_conditions/cbvco2", 0.);
	stateVars[26] = ifile("initial_conditions/cvo2", 0.);
	stateVars[27] = ifile("initial_conditions/cvco2", 0.);
	}

	HsysArt = ifile("initial_conditions/H", 0.);
	HpulCap = ifile("initial_conditions/H", 0.);
	pHsysArt = ifile("initial_conditions/pH", 0.);
	pHpulCap = ifile("initial_conditions/pH", 0.);
	AlveoliLungsAlg[7] = 90.0;
	AlveoliLungsAlg[8] = 40.0;
	double cCO2_old =  0.535371/22.7; // 22.0e-3; // mol/l

	if (initial_cond = 1) {
		// leggi le pressioni poi con dissociation calcolo le concentrazioni e inizializzale
		// prova se con _old funziona tutto bene, poi prova con dissociation nuovo, inizializza 
		// anche po2 co2 sistemica, vedendo quanto vale la concentrazione cao2, caco2 (var alg)
		PressuresInitialCond[0] = ifile("initial_conditions/pppco2", 38.9974);
		PressuresInitialCond[1] = ifile("initial_conditions/pppo2", 104.082);
		PressuresInitialCond[2] = ifile("initial_conditions/pepco2", 37.6065);
		PressuresInitialCond[3] = ifile("initial_conditions/pepo2", 67.5991);
		PressuresInitialCond[4] = ifile("initial_conditions/pspco2", 45.6982);
		PressuresInitialCond[5] = ifile("initial_conditions/pspo2", 30.6311);
		PressuresInitialCond[6] = ifile("initial_conditions/pmpco2", 40.7813);
		PressuresInitialCond[7] = ifile("initial_conditions/pmpo2", 38.0859);
		PressuresInitialCond[8] = ifile("initial_conditions/phpco2", 49.6231);
		PressuresInitialCond[9] = ifile("initial_conditions/phpo2", 40.5419);
		PressuresInitialCond[10] = ifile("initial_conditions/pbpco2", 43.9942);
		PressuresInitialCond[11] = ifile("initial_conditions/pbpo2", 32.7514);
		PressuresInitialCond[12] = ifile("initial_conditions/pevco2", 37.6065);
		PressuresInitialCond[13] = ifile("initial_conditions/pevo2", 67.5991);
		PressuresInitialCond[14] = ifile("initial_conditions/psvco2", 45.6982);
		PressuresInitialCond[15] = ifile("initial_conditions/psvo2", 30.6311);
		PressuresInitialCond[16] = ifile("initial_conditions/pmvco2", 40.7813);
		PressuresInitialCond[17] = ifile("initial_conditions/pmvo2", 38.0859);
		PressuresInitialCond[18] = ifile("initial_conditions/phvco2", 49.6231);
		PressuresInitialCond[19] = ifile("initial_conditions/phvo2", 40.5419);
		PressuresInitialCond[20] = ifile("initial_conditions/pbvco2", 43.9942);
		PressuresInitialCond[21] = ifile("initial_conditions/pbvo2", 32.7514);
		PressuresInitialCond[22] = ifile("initial_conditions/pvco2", 42.9239);
		PressuresInitialCond[23] = ifile("initial_conditions/pvo2", 39.0809);

		AlveoliLungsAlg[3] = PressuresInitialCond[1];
		AlveoliLungsAlg[4] = PressuresInitialCond[0];

		for (int i = 0; i < 12; ++i) {
			double pCO2 = PressuresInitialCond[2 * i] * 0.133322;
			double pO2 = PressuresInitialCond[2 * i + 1] * 0.133322;
			tie(stateVars[2 * i + 4], stateVars[2 * i + 5]) = dissociation(pCO2, pO2, cCO2_old, HsysArt);

			stateVars[2 * i + 4] *= 22.7;
			stateVars[2 * i + 5] *= 22.7;
		}

		// cout << "Initial conditions from pressures: " << endl;
		// cout << "pppco2 = " << PressuresInitialCond[0] << ", pppo2 = " << PressuresInitialCond[1] << endl;

		// cout << "Initial conditions from pressures: " << endl;
		// cout << "cppo2 = " << stateVars[4] << ", cppco2 = " << stateVars[5] << endl;
		// cout << "cepo2 = " << stateVars[6] << ", cepco2 = " << stateVars[7] << endl;
		// cout << "cspo2 = " << stateVars[8] << ", cspco2 = " << stateVars[9] << endl;
		// cout << "cmpo2 = " << stateVars[10] << ", cmpco2 = " << stateVars[11] << endl;
		// cout << "chpo2 = " << stateVars[12] << ", chpco2 = " << stateVars[13] << endl;
		// cout << "cbpo2 = " << stateVars[14] << ", cbpco2 = " << stateVars[15] << endl;
		// cout << "cevo2 = " << stateVars[16] << ", cevco2 = " << stateVars[17] << endl;
		// cout << "csvo2 = " << stateVars[18] << ", csvco2 = " << stateVars[19] << endl;
		// cout << "cmvo2 = " << stateVars[20] << ", cmvco2 = " << stateVars[21] << endl;
		// cout << "chvo2 = " << stateVars[22] << ", chvco2 = " << stateVars[23] << endl;
		// cout << "cbvo2 = " << stateVars[24] << ", cbvco2 = " << stateVars[25] << endl;
		// cout << "cvo2 = " << stateVars[26] << ", cvco2 = " << stateVars[27] << endl;

		// tie(stateVars[4], stateVars[5]) = dissociation_old(PressuresInitialCond[1], PressuresInitialCond[0]); 
		// tie(stateVars[6],  stateVars[7])  = dissociation_old(PressuresInitialCond[3],  PressuresInitialCond[2]); 
		// tie(stateVars[8],  stateVars[9])  = dissociation_old(PressuresInitialCond[5],  PressuresInitialCond[4]); 
		// tie(stateVars[10], stateVars[11]) = dissociation_old(PressuresInitialCond[7],  PressuresInitialCond[6]); 
		// tie(stateVars[12], stateVars[13]) = dissociation_old(PressuresInitialCond[9],  PressuresInitialCond[8]); 
		// tie(stateVars[14], stateVars[15]) = dissociation_old(PressuresInitialCond[11], PressuresInitialCond[10]); 
		// tie(stateVars[16], stateVars[17]) = dissociation_old(PressuresInitialCond[13], PressuresInitialCond[12]); 
		// tie(stateVars[18], stateVars[19]) = dissociation_old(PressuresInitialCond[15], PressuresInitialCond[14]);
		// tie(stateVars[20], stateVars[21]) = dissociation_old(PressuresInitialCond[17], PressuresInitialCond[16]);
		// tie(stateVars[22], stateVars[23]) = dissociation_old(PressuresInitialCond[19], PressuresInitialCond[18]);
		// tie(stateVars[24], stateVars[25]) = dissociation_old(PressuresInitialCond[21], PressuresInitialCond[20]);
		// tie(stateVars[26], stateVars[27]) = dissociation_old(PressuresInitialCond[23], PressuresInitialCond[22]);
	}
	if (verbose){
		cout << "Printing initial conditions State Vars: ... " << endl;
		for (int i = 0; i < stateVars.size(); i++) {    
			cout << "State Vars" << i << " : " << stateVars[i] << endl;     
		}
	}
	
	getAlgebraicRelations();
}

void gasTransport::printGASparameters() {
	cout << "[Lung Gas Parama]" << endl;
	for (int i = 0; i < LungsGasParams.size(); i++) {
		cout << i << " value: " << LungsGasParams[i] << endl;
	}
	cout << "[Tissues Venous Gas Params]" << endl;
	for (int i = 0; i < TissuesVenousGasParam.size(); i++) {
		cout << i << " value: " << TissuesVenousGasParam[i] << endl;
	}
	cout << "delayed params: " << DelayedGasParam[0] << " " << DelayedGasParam[1] << endl;
}


tuple<double, double> gasTransport::calculate_pH(double cCO2, double H_old) {
    /**
     * Risolve l’equazione cubica g([H+]) = 0 con Newton–Raphson e restituisce pH = -log10([H+]) 
     */

    // Costanti fisiologiche
    const double KaCO2 = pow(10, -6.1);      // [mol/l] costante dissociazione CO2 (42)
    const double KaPr = pow(10, -7.3);       // [mol/l] costante dissociazione Hb (43)
    const double NaOH_0 = 46.2e-3;     		  // [mol/l] [NaOH] iniziale (46)
    const double HPr_0 = 39.8e-3;      		  // [mol/l] [HPr] iniziale (47)

	const double tol = 1e-4 * H_old;     	  // tolleranza
	const int max_iter = 100;

    // Coefficienti eq. cubica
    double a2 = KaPr + NaOH_0 + KaCO2;
    double a1 = KaCO2 * (NaOH_0 - cCO2) + KaPr * (KaCO2 + NaOH_0 - HPr_0);
    double a0 = KaCO2 * KaPr * (NaOH_0 - HPr_0 - cCO2);

    // Funzione g(H+) e sua derivata
    auto g = [&](double H) {
        return pow(H, 3) + a2 * pow(H, 2) + a1 * H + a0;  // (48, 61)
    };
    auto dg = [&](double H) {
        return 3 * pow(H, 2) + 2 * a2 * H + a1;
    };

    // Newton–Raphson (62)
    double H = H_old;  // initial value [H+] mol/L
    for (int i = 0; i < max_iter; ++i) {
        double delta = g(H) / dg(H);
        H -= delta;
        if (fabs(delta) < tol) break;
    }

    return make_tuple(-log10(H), H);  // pH
}

double gasTransport::calculate_SO2(double pO2, double pCO2, double pH) {
	/**
     * Saturation of O2 in blood
     */
	double T_C = 37.0; 		// [°C] blood temperature in Celsius
	double xHbf = 0.0; 		// [-] fraction of fetal hemoglobin concentration in the blood
	double xHbCO = 0.005; 	// [-] fraction of carboxyhemoglobin in the blood
	double xHi = 0.005; 	// [-] fraction of glycohemoglobin in the blood
	double cDPG = 5.01e-3 ; // [mol/l] concentration of 2,3-diphosphogycerate in the erythrocyte

    double a = -0.72 * (pH - 7.4)										 // (27)
               + 0.09 * log(pCO2 / 5.33)
               + (0.07 - 0.03 * xHbf) * (cDPG - 5.0e-3)
               - 0.368 * xHbCO
               - 0.174 * xHi
               - 0.28 * xHbf;

    double x0 = 1.946 + a + 0.055 * (T_C - 37.0); 					 // (26)
    double x = log(pO2); 									 // (25)
    double h = 3.5 + a; 											 // (24)	
    double y = 1.875 + x - x0 + h * tanh(0.5343 * (x - x0));		 // (23)

    double SO2 = 1.0 / (1.0 + exp(-y)); 							 // (22)
	return SO2;
}


double gasTransport::calculate_cO2(double pO2, double SO2) {

	static bool first_call = true;
    /**
     * Total O2 concentration in blood = free + concentration in red blood cells
     */
	double alphaO2 = 9.83e-6;  	// [mol/l/kPa] solubility of O2
	double cHb = 9.3e-3; 		// [mol/l] total concentration of Hb 

	double cO2 = alphaO2 * pO2 + cHb * SO2; 		// (21) 
    return cO2; 									// cO2 esce in mol/l
}		 
	
double gasTransport::calculate_cCO2(double pCO2, double pH, double SO2) {
    /**
     *Concentration of CO2 in blood (eritrociti + plasma)  
     */
	double alphaCO2_Ery = 0.195e-3;  	// [mol/l/kPA] solubility of CO2 in erythrocyte (33)
	double alphaCO2_Pla = 0.230e-3; 	// [mol/l/kPA] solubility of CO2 in plasma (34)
	double cHb = 9.3e-3; 				// [mol/l] total concentration of Hb 
	double cHb_Ery = 21.0e-3; 			// [mol/l] concentration of Hb in erythrocyte
	double cHbratio = cHb/cHb_Ery; 		// hemoglobin concentration of the blood (29)
    
	// pK erythrocyte e plasma
    double pHEry = 7.19 + 0.77 * (pH - 7.4) + 0.035 * (1.0 - SO2);					// (36)
	double pKEry = 6.125 - log10(1.0 + pow(10.0, pHEry - 7.84 - 0.06 * SO2)); 		// (35)
    double pKPla = 6.125 - log10(1.0 + pow(10.0, pH - 8.7)); //+ alphaCO2_Pla;						// (37)

    // cCO2 erythrocyte e plasma
    double cCO2_Ery = alphaCO2_Ery * pCO2 * (1.0 + pow(10.0, pHEry - pKEry));		// (33)
    double cCO2_Pla = alphaCO2_Pla * pCO2 * (1.0 + pow(10.0, pH - pKPla));			// (34)

    double cCO2 = (cCO2_Ery * cHbratio + cCO2_Pla * (1.0 - cHbratio)); 				// (29)
    return cCO2;						// cO2 esce in mol/l
}

double gasTransport::dpH_dcCO2(double cCO2, double H_old) {
	double pH = std::get<0>(calculate_pH(cCO2, H_old));
    const double h = 1e-4 * pH;
    double pH1 = std::get<0>(calculate_pH(cCO2 + h, H_old));
    double pH0 = std::get<0>(calculate_pH(cCO2 - h, H_old));
    return (pH1 - pH0) / (2 * h);
}

double gasTransport::dphiCO2_dpH(double pCO2, double pH, double SO2) {
	double phi = calculate_cCO2(pCO2, pH, SO2);
    const double h = 1e-4 * phi;
    double phi1 = calculate_cCO2(pCO2, pH + h, SO2);
    double phi0 = calculate_cCO2(pCO2, pH - h, SO2);
    return (phi1 - phi0) / (2 * h);
}

double gasTransport::df_dcCO2(double cCO2, double pCO2, double SO2, double H_old) {
    double pH = std::get<0>(calculate_pH(cCO2, H_old));
    double dphi_dpH = dphiCO2_dpH(pCO2, pH, SO2);
    double dpH_dc = dpH_dcCO2(cCO2, H_old);
    return dphi_dpH * dpH_dc - 1.0;
}

tuple<double, double> gasTransport::dissociation(double pCO2, double pO2, double cCO2_old, double H_old) { 
    const double tol = 1e-4 * cCO2_old;
    const int max_iter = 100;

    double cCO2 = cCO2_old; // 22.0e-3 mol/l = 0.4994 mLgas/mLblood;  // iniziale

    for (int i = 0; i < max_iter; ++i) {
        double pH = get<0>(calculate_pH(cCO2, H_old));
		double SO2 = calculate_SO2(pO2, pCO2, pH);
        double phi = calculate_cCO2(pCO2, pH, SO2);
        double f = phi - cCO2;
        double df = df_dcCO2(cCO2, pCO2, SO2, H_old);
		
		// Newton-Raphson step
		double delta = f / df;
        cCO2 -= delta;

        if (fabs(delta) < tol) break;
    }

    double pH_out = get<0>(calculate_pH(cCO2, H_old));
	double SO2_out = calculate_SO2(pO2, pCO2, pH_out);
    double cO2_out = calculate_cO2(pO2, SO2_out); 
    double cCO2_out = cCO2;

	return make_tuple(cO2_out, cCO2_out);
}


tuple<double, double> gasTransport::dissociation_old(double p1, double p2){
	/**
	* takes PO2(p1) and PCO2(p2) partial pressures and provides concentration of O2 and CO2 in blood,
	* as in Spencer(1979), eqs(4) and (5)
	*/
	double f1;
	double f2;
	double c1;
	double c2;
	double LungsGasParams_5= 0.20448;		// csat1
	double LungsGasParams_6= 1.9564192;		// csat2
	double LungsGasParams_7= 0.3836;		// h1
	double LungsGasParams_8= 1.819;			// h2
	double LungsGasParams_9= 0.03198; 		// alpha1
	double LungsGasParams_10= 0.05591; 		// alpha2
	double LungsGasParams_11 = 0.008275; 	// beta1
	double LungsGasParams_12 = 0.03255;		// beta2
	double LungsGasParams_13= 14.99; 		// k1
	double LungsGasParams_14= 194.4; 		// k2


	f1 = p1 * (1. + LungsGasParams_11 * p2) / LungsGasParams_13 / (1. + LungsGasParams_9 * p2);
	f2 = p2 * (1. + LungsGasParams_12 * p1) / LungsGasParams_14 / (1. + LungsGasParams_10 * p1);

	c1 = LungsGasParams_5 * pow(f1,(1. / LungsGasParams_7)) / (1. + pow(f1,(1. / LungsGasParams_7)));
	c2 = LungsGasParams_6 * pow(f2,(1. / LungsGasParams_8)) / (1. + pow(f2,(1. / LungsGasParams_8)));

	return make_tuple(c1, c2);
}

tuple<double, double> gasTransport::invertedDissociation_old(double c1, double c2){
	/**
	* takes O2 and CO2 concentration and provides partial O2 and CO2 blood pressures,
	* as in Spencer(1979), eqs(6) and (7)
	*/
	double d1;
	double d2;
	double r1;
	double r2;
	double s1;
	double s2;
	double p1;
	double p2;
	double LungsGasParams_5= 0.20448;		// csat1
	double LungsGasParams_6= 1.9564192;		// csat2
	double LungsGasParams_7= 0.3836;		// h1
	double LungsGasParams_8= 1.819;			// h2
	double LungsGasParams_9= 0.03198; 		// alpha1
	double LungsGasParams_10= 0.05591; 		// alpha2
	double LungsGasParams_11 = 0.008275; 	// beta1
	double LungsGasParams_12 = 0.03255;		// beta2
	double LungsGasParams_13= 14.99; 		// k1
	double LungsGasParams_14= 194.4; 		// k2

	d1 = max(LungsGasParams_13*pow((c1 / (LungsGasParams_5 - c1)),LungsGasParams_7), 0.);
	d2 = max(LungsGasParams_14*pow((c2 / (LungsGasParams_6 - c2)),LungsGasParams_8), 0.);

	r1 = -1.*(1. + LungsGasParams_11 * d2 - LungsGasParams_12 * d1 - LungsGasParams_9 * LungsGasParams_10*d1*d2) / 2. / (LungsGasParams_12 + LungsGasParams_10 * LungsGasParams_11*d2);
	r2 = -1.*(1. + LungsGasParams_12 * d1 - LungsGasParams_11 * d2 - LungsGasParams_10 * LungsGasParams_9*d2*d1) / 2. / (LungsGasParams_11 + LungsGasParams_9 * LungsGasParams_12*d1);

	s1 = -1.*(d1 + LungsGasParams_9 * d1*d2) / (LungsGasParams_12 + LungsGasParams_10 * LungsGasParams_11*d2);
	s2 = -1.*(d2 + LungsGasParams_10 * d2*d1) / (LungsGasParams_11 + LungsGasParams_9 * LungsGasParams_12*d1);

	p1 = r1 + pow((pow(r1,2) - s1),0.5);
	p2 = d2 * (1 + LungsGasParams_10 * p1) / (1 + LungsGasParams_12 * p1);
	
	return make_tuple(p1, p2);
}



// tuple<double, double> gasTransport::invertedDissociation(double cO2_target_mL, double cCO2_target_mL, double H_old, double pO2_OLD, double pCO2_OLD) {

//     double pCO2_mmHg = pCO2_OLD;   // stima iniziale
//     double pO2_mmHg = pO2_OLD;     // stima iniziale

// 	// convertire concentrazioni da mLgas/mLblood a mol/l
// 	double cCO2_target = cCO2_target_mL / 22.7; 
// 	double cO2_target = cO2_target_mL / 22.7;
// 	cout << "cCO2_target = " << cCO2_target << ", cO2_target = " << cO2_target << endl;
// 	double cCO2_old = cCO2_target; // stima iniziale
// 	double cO2_old = cO2_target; // stima iniziale
// 	// convertire pressioni da mmHg a kPa 
// 	double pO2 =  pO2_mmHg * 0.133322; 
// 	double pCO2 =  pCO2_mmHg * 0.133322;

// 	const double hCO2 = 1e-4 * pCO2;
// 	const double hO2 = 1e-2 * pO2;
//     const double tolCO2 = 1e-3 * pCO2;
// 	const double tolO2 = 1e-3 * pO2;
//     const int max_iter = 100;
// 	double cO2, cCO2;
// 	double cO2_pluspO2, cCO2_pluspO2;
// 	double cO2_minuspO2, cCO2_minuspO2;
// 	double cO2_pluspCO2, cCO2_pluspCO2;
// 	double cO2_minuspCO2, cCO2_minuspCO2;


//     for (int iter = 0; iter < max_iter; ++iter) {
//     	tie(cO2, cCO2) = dissociation(pCO2, pO2, cCO2_old, H_old);
// 		cout << "cO2 = " << cO2 << ", cCO2 = " << cCO2 << endl;
// 		tie(cO2_pluspO2, cCO2_pluspO2) = dissociation(pCO2, pO2 + hO2, cCO2_old, H_old);
// 		tie(cO2_minuspO2, cCO2_minuspO2) = dissociation(pCO2, pO2 - hO2, cCO2_old, H_old);
// 		tie(cO2_pluspCO2, cCO2_pluspCO2) = dissociation(pCO2 + hCO2, pO2, cCO2_old, H_old);
// 		tie(cO2_minuspCO2, cCO2_minuspCO2) = dissociation(pCO2 - hCO2, pO2, cCO2_old, H_old);

// 		double f1 = cO2 - cO2_target;
//         double f2 = cCO2 - cCO2_target;

//         // Derivate parziali numeriche
//         double df1_pO2   = (cO2_pluspO2 - cO2_minuspO2) / (2*hO2);
//         double df1_pCO2  = (cO2_pluspCO2 - cO2_minuspCO2) / (2*hCO2);
//         double df2_pO2  = (cCO2_pluspO2 - cCO2_minuspO2) / (2*hO2);
//         double df2_pCO2 = (cCO2_pluspCO2 - cCO2_minuspCO2) / (2*hCO2);
// 		cout << "f1 = " << f1 << ", f2 = " << f2 << endl;
// 		cout << "df1_pO2 = " << df1_pO2 << ", df1_pCO2 = " << df1_pCO2 << endl;
// 		cout << "df2_pO2 = " << df2_pO2 << ", df2_pCO2 = " << df2_pCO2 << endl;

//         // Matrice Jacobiana
//         double det = df2_pCO2 * df1_pO2 - df2_pO2 * df1_pCO2;
//         if (fabs(det) < 1e-5){
// 			cout<<"Determinant is too small, division by zero avoided."<<endl;
// 		}; 
// 		cout << "det = " << det << endl;

//         // Inversa della Jacobiana * [f1, f2]
//         double dpO2  = (f1 * df2_pCO2 / det - f2 * df1_pCO2 / det);
//         double dpCO2 = (-f1 * df2_pO2 / det + f2 * df1_pO2 / det);

// 		cout << "dpO2 = " << dpO2 << ", dpCO2 = " << dpCO2 << endl;

//         pCO2 -= dpCO2;
//         pO2  -= dpO2;

//         if (fabs(dpCO2) < tolCO2 && fabs(dpO2) < tolO2)
//             break;
//     }

// 	// convertire pressioni da kPa a mmHg
// 	double pO2_out =  pO2 / 0.133322; 
// 	double pCO2_out =  pCO2 / 0.133322;
// 	cout << "pO2_out = " << pO2_out << ", pCO2_out = " << pCO2_out << endl;

//     return make_tuple(pCO2_out, pO2_out);
// }


// tuple<double, double> gasTransport::invertedDissociation(double cO2_target_mL, double cCO2_target_mL, double H_old, double pO2_OLD, double pCO2_OLD) {

//     double pCO2_mmHg = pCO2_OLD;   // stima iniziale
//     double pO2_mmHg = pO2_OLD;     // stima iniziale

// 	// convertire concentrazioni da mLgas/mLblood a mol/l
// 	double cCO2_target = cCO2_target_mL / 22.7; 
// 	double cO2_target = cO2_target_mL / 22.7;
// 	//cout << "cCO2_target = " << cCO2_target << ", cO2_target = " << cO2_target << endl;
// 	double cCO2_old = cCO2_target; // stima iniziale
// 	double cO2_old = cO2_target; // stima iniziale
// 	// convertire pressioni da mmHg a kPa 
// 	double pO2 =  pO2_mmHg * 0.133322; 
// 	double pCO2 =  pCO2_mmHg * 0.133322;

// 	const double hCO2 = 1e-4 * pCO2;
// 	const double hO2 = 1e-2 * pO2;
//     const double tolCO2 = 1e-3 * pCO2;
// 	const double tolO2 = 1e-3 * pO2;
//     const int max_iter = 100;
// 	double cO2, cCO2;
// 	double cO2_pluspO2, cCO2_pluspO2;
// 	double cO2_minuspO2, cCO2_minuspO2;
// 	double cO2_pluspCO2, cCO2_pluspCO2;
// 	double cO2_minuspCO2, cCO2_minuspCO2;

// 	// Parametro di regolarizzazione
// 	const double lambda = 1e-5; 

//     for (int iter = 0; iter < max_iter; ++iter) {
//     	tie(cO2, cCO2) = dissociation(pCO2, pO2, cCO2_old, H_old);
// 		// cout << "cO2 = " << cO2 << ", cCO2 = " << cCO2 << endl;
// 		tie(cO2_pluspO2, cCO2_pluspO2) = dissociation(pCO2, pO2 + hO2, cCO2_old, H_old);
// 		tie(cO2_minuspO2, cCO2_minuspO2) = dissociation(pCO2, pO2 - hO2, cCO2_old, H_old);

// 		tie(cO2_pluspCO2, cCO2_pluspCO2) = dissociation(pCO2 + hCO2, pO2, cCO2_old, H_old);
// 		tie(cO2_minuspCO2, cCO2_minuspCO2) = dissociation(pCO2 - hCO2, pO2, cCO2_old, H_old);

// 		double f1 = cO2 - cO2_target;
//         double f2 = cCO2 - cCO2_target;

//         // Derivate parziali numeriche
//         double df1_pO2   = (cO2_pluspO2 - cO2_minuspO2) / (2*hO2);
//         double df1_pCO2  = (cO2_pluspCO2 - cO2_minuspCO2) / (2*hCO2);
//         double df2_pO2  = (cCO2_pluspO2 - cCO2_minuspO2) / (2*hO2);
//         double df2_pCO2 = (cCO2_pluspCO2 - cCO2_minuspCO2) / (2*hCO2);
// 		// cout << "f1 = " << f1 << ", f2 = " << f2 << endl;
// 		// cout << "df1_pO2 = " << df1_pO2 << ", df1_pCO2 = " << df1_pCO2 << endl;
// 		// cout << "df2_pO2 = " << df2_pO2 << ", df2_pCO2 = " << df2_pCO2 << endl;

//         // Matrice Jacobiana con regolarizzazione
//         double det = df2_pCO2 * df1_pO2 - df2_pO2 * df1_pCO2;
        
//         // Aggiungi un termine di regolarizzazione alla matrice Jacobiana
//         double regularized_det = det + lambda;

//         // Verifica che il determinante regolarizzato non sia troppo vicino a zero
//         if (fabs(regularized_det) < 1e-5){
// 			// cout << "Determinant is too small, division by zero avoided." << endl;
// 		}; 
// 		// cout << "regularized_det = " << regularized_det << endl;

//         // Inversa della Jacobiana regolarizzata * [f1, f2]
//         double dpO2  = (f1 * df2_pCO2 / regularized_det - f2 * df1_pCO2 / regularized_det);
//         double dpCO2 = (-f1 * df2_pO2 / regularized_det + f2 * df1_pO2 / regularized_det);

// 		// cout << "dpO2 = " << dpO2 << ", dpCO2 = " << dpCO2 << endl;

//         pCO2 -= dpCO2;
//         pO2  -= dpO2;

//         if (fabs(dpCO2) < tolCO2 && fabs(dpO2) < tolO2)
// 			//cout << "Converged after " << iter << " iterations." << endl;
//             break;
//     }

// 	// convertire pressioni da kPa a mmHg
// 	double pO2_out =  pO2 / 0.133322; 
// 	double pCO2_out =  pCO2 / 0.133322;
// 	// cout << "pO2_out = " << pO2_out << ", pCO2_out = " << pCO2_out << endl;

//     return make_tuple(pCO2_out, pO2_out);
// }

tuple<double, double> gasTransport::invertedDissociation(double cO2_target_mL, double cCO2_target_mL, double H_old, double pO2_OLD, double pCO2_OLD) {

    double pCO2_mmHg = pCO2_OLD;   // stima iniziale
    double pO2_mmHg = pO2_OLD;     // stima iniziale

	// convertire concentrazioni da mLgas/mLblood a mol/l
	double cCO2_target = cCO2_target_mL / 22.7; 
	double cO2_target = cO2_target_mL / 22.7;
	// cout << "cCO2_target = " << cCO2_target << ", cO2_target = " << cO2_target << endl;
	double cCO2_old = cCO2_target; // stima iniziale
	double cO2_old = cO2_target; // stima iniziale
	// convertire pressioni da mmHg a kPa 
	double pO2 =  pO2_mmHg * 0.133322; 
	double pCO2 =  pCO2_mmHg * 0.133322;

	const double hCO2 = 1e-4 * pCO2;
	const double hO2 = 1e-4 * pO2;
    const double tol = 1e-6;
    const int max_iter = 100;
	double cO2, cCO2;
	double cO2_pluspO2, cCO2_pluspO2;
	double cO2_minuspO2, cCO2_minuspO2;
	double cO2_pluspCO2, cCO2_pluspCO2;
	double cO2_minuspCO2, cCO2_minuspCO2;

    for (int iter = 0; iter < max_iter; ++iter) {
		// cout << "iter = " << iter << endl;
    	tie(cO2, cCO2) = dissociation(pCO2, pO2, cCO2_old, H_old);
		// cout << "cO2 = " << cO2 << ", cCO2 = " << cCO2 << endl;
		tie(cO2_pluspO2, cCO2_pluspO2) = dissociation(pCO2, pO2 + hO2, cCO2_old, H_old);
		tie(cO2_minuspO2, cCO2_minuspO2) = dissociation(pCO2, pO2 - hO2, cCO2_old, H_old);
		tie(cO2_pluspCO2, cCO2_pluspCO2) = dissociation(pCO2 + hCO2, pO2, cCO2_old, H_old);
		tie(cO2_minuspCO2, cCO2_minuspCO2) = dissociation(pCO2 - hCO2, pO2, cCO2_old, H_old);

		double f1 = cO2 - cO2_target;
        double f2 = cCO2 - cCO2_target;
		// cout << cCO2_target << " " << cCO2 << endl;

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

		// cout << "dpO2 = " << dpO2 << ", dpCO2 = " << dpCO2 << endl;

        double alpha = 1.0;
		// for (int j = 0; j < max_iter; ++j) {
		// 	double new_pO2 = pO2 - alpha * dpO2;
		// 	double new_pCO2 = pCO2 - alpha * dpCO2;
		// 	double new_cO2, new_cCO2;

		// 	tie(new_cO2, new_cCO2) = dissociation(new_pCO2, new_pO2, cCO2_old, H_old);
		// 	double f1_new = new_cO2 - cO2_target;
		// 	double f2_new = new_cCO2 - cCO2_target;

		// 	double norm_f_new = std::max(std::abs(f1_new), std::abs(f2_new));
		// 	double norm_f_old = std::max(std::abs(f1), std::abs(f2));

		// 	if (norm_f_new > norm_f_old) break; // range pressioni
			
		// 	alpha *= 0.8;
		// }

		pCO2 -= alpha * dpCO2;
		pO2  -= alpha * dpO2;

        double rel_err = std::max(abs(f1) / cO2_target, abs(f2) / cCO2_target);
		if (rel_err < tol)  break;
    }

	// convertire pressioni da kPa a mmHg
	double pO2_out =  pO2 / 0.133322; 
	double pCO2_out =  pCO2 / 0.133322;
	// cout << "pO2_out = " << pO2_out << ", pCO2_out = " << pCO2_out << endl;

    return make_tuple(pCO2_out, pO2_out);
}

double gasTransport::heaviside(double value1, double value2){
	double value;
	if (value1 < 0) {
		value = 0;
	}
	if (value1 == 0) {
		value = value2;
	}
	if (value1 > 0) {
		value = 1;
	}
	return value;
}

void gasTransport::getAlgebraicRelations() {
	AlveoliLungsAlg[0] = heaviside(CommonVars[1], 0.);      //heaviside from vdot
	AlveoliLungsAlg[1] = heaviside(-1.*CommonVars[1], 0.);
		
	AlveoliLungsAlg[12] = max(stateVars[2], 0.)*(LungsGasParams[3] - LungsGasParams[4]);  // A - pAo2
	AlveoliLungsAlg[2] = max(stateVars[3], 0.)*(LungsGasParams[3] - LungsGasParams[4]); // A53 - pAco2

	tie (pHpulCap, HpulCap) = calculate_pH(stateVars[5]/ 22.7, HpulCap); // pH	
	// cout << "Here0" << endl;
	// cout << stateVars[4] << " " << stateVars[5] << " " << HpulCap << " " << pHpulCap << " " << AlveoliLungsAlg[3] << " " << AlveoliLungsAlg[4] << endl;
	tie(AlveoliLungsAlg[4], AlveoliLungsAlg[3]) = invertedDissociation(stateVars[4], stateVars[5], HpulCap, AlveoliLungsAlg[3], AlveoliLungsAlg[4]); // pppo2,pppco2
	// cout << AlveoliLungsAlg[4] << " " << AlveoliLungsAlg[3] << endl;
	double denom = CommonVars[31] + CommonVars[32];
	if (denom != 0) {
		AlveoliLungsAlg[5] = max((CommonVars[31]*stateVars[4] + CommonVars[32]*TissueVenousGasAlg[2]) / denom, 0.0); //A54 - cao2
		AlveoliLungsAlg[6] = max((CommonVars[31]*stateVars[5] + CommonVars[32]*TissueVenousGasAlg[3]) / denom, 0.0); //A55 - caco2
	} else {
		AlveoliLungsAlg[5] = 0.212045;  // mLgas/mLblood
		AlveoliLungsAlg[6] = 0.535371;  // mLgas/mLblood
	}

	tie (pHsysArt, HsysArt) = calculate_pH(AlveoliLungsAlg[6]/ 22.7, HsysArt);
	// cout << "Here1" << endl;
	// cout << AlveoliLungsAlg[5] << " " << AlveoliLungsAlg[6] << " " << HsysArt << " " << pHsysArt << " " << AlveoliLungsAlg[7] << " " << AlveoliLungsAlg[8] << endl;
	tie(AlveoliLungsAlg[8], AlveoliLungsAlg[7]) = invertedDissociation(AlveoliLungsAlg[5], AlveoliLungsAlg[6], HsysArt, AlveoliLungsAlg[7], AlveoliLungsAlg[8]); //pao2,paco2
	// cout << AlveoliLungsAlg[8] << " " << AlveoliLungsAlg[7] << endl;
	if ((LungsGasParams[5] - LungsGasParams[16]*1.34 - AlveoliLungsAlg[7]*0.003 / 100.) >= 0.) {
		AlveoliLungsAlg[9] = AlveoliLungsAlg[7] * 0.003 / 100.; //cao2diss
	}
	else {
		AlveoliLungsAlg[9] = LungsGasParams[5] - LungsGasParams[16] * 1.34;  //cao2diss
	}
	AlveoliLungsAlg[10] = LungsGasParams[17] * (AlveoliLungsAlg[12] - AlveoliLungsAlg[3]);  //mdoto2
	AlveoliLungsAlg[11] = LungsGasParams[18] * (AlveoliLungsAlg[2] - AlveoliLungsAlg[4]);  //mdotco2
}

void gasTransport::getTimeDerivative() {
	getAlgebraicRelations();
	// lung exchange
	dvdt[0] = 1. / CommonVars[0]*(AlveoliLungsAlg[0] * CommonVars[1]*(LungsGasParams[0] - max(stateVars[0], 0.)) + AlveoliLungsAlg[1] * CommonVars[2]*(max(stateVars[0], 0.) - max(stateVars[2], 0.))); // A42
	dvdt[1] = 1. / CommonVars[0]*(AlveoliLungsAlg[0] * CommonVars[1]*(LungsGasParams[1] - max(stateVars[1], 0.)) + AlveoliLungsAlg[1] * CommonVars[2]*(max(stateVars[1], 0.) - max(stateVars[3], 0.))); // A43
	dvdt[2] = 1. / CommonVars[3]*(AlveoliLungsAlg[0] * CommonVars[2]*(max(stateVars[0], 0.) - max(stateVars[2], 0.)) - LungsGasParams[2] * (AlveoliLungsAlg[10])); // A44
	dvdt[3] = 1. / CommonVars[3]*(AlveoliLungsAlg[0] * CommonVars[2]*(max(stateVars[1], 0.) - max(stateVars[3], 0.)) - LungsGasParams[2] * (AlveoliLungsAlg[11]));  // A45
	dvdt[4] = 1. / CommonVars[13]*(CommonVars[31]*(TissueVenousGasAlg[2] - stateVars[4]) + AlveoliLungsAlg[10]);
	dvdt[5] = 1. / CommonVars[13]*(CommonVars[31]*(TissueVenousGasAlg[3] - stateVars[5]) + AlveoliLungsAlg[11]);

	// tissue gas exchange
	dvdt[6] = 1. / (TissuesVenousGasParam[0] + CommonVars[4])*(CommonVars[26]*(TissueVenousGasAlg[0] - max(stateVars[6], 0.)) - TissuesVenousGasParam[1]); // A63
	dvdt[7] = 1. / (TissuesVenousGasParam[0] + CommonVars[4])*(CommonVars[26]*(TissueVenousGasAlg[1] - max(stateVars[7], 0.)) + TissuesVenousGasParam[2]); // A64
	dvdt[8] = 1. / (TissuesVenousGasParam[3] + CommonVars[5])*(CommonVars[27]*(TissueVenousGasAlg[0] - max(stateVars[8], 0.)) - TissuesVenousGasParam[4]); // A65
	dvdt[9] = 1. / (TissuesVenousGasParam[3] + CommonVars[5])*(CommonVars[27]*(TissueVenousGasAlg[1] - max(stateVars[9], 0.)) + TissuesVenousGasParam[5]);  // A66
	dvdt[10] = 1. / (TissuesVenousGasParam[6] + CommonVars[6])*(CommonVars[28]*(TissueVenousGasAlg[0] - max(stateVars[10], 0.)) - TissuesVenousGasParam[7]); // A61
	dvdt[11] = 1. / (TissuesVenousGasParam[6] + CommonVars[6])*(CommonVars[28]*(TissueVenousGasAlg[1] - max(stateVars[11], 0.)) + TissuesVenousGasParam[8]); // A62
	dvdt[12] = 1. / (TissuesVenousGasParam[9] + CommonVars[7])*(CommonVars[29]*(TissueVenousGasAlg[0] - max(stateVars[12], 0.)) - TissuesVenousGasParam[10]); // A57
	dvdt[13] = 1. / (TissuesVenousGasParam[9] + CommonVars[7])*(CommonVars[29]*(TissueVenousGasAlg[1] - max(stateVars[13], 0.)) + TissuesVenousGasParam[11]); // A58
	dvdt[14] = 1. / (TissuesVenousGasParam[12] + CommonVars[8])*(CommonVars[30]*(TissueVenousGasAlg[0] - max(stateVars[14], 0.)) - TissuesVenousGasParam[13]); // A59
	dvdt[15] = 1. / (TissuesVenousGasParam[12] + CommonVars[8])*(CommonVars[30]*(TissueVenousGasAlg[1] - max(stateVars[15], 0.)) + TissuesVenousGasParam[14]); // A60

	// venous pool gas transport
	dvdt[16] = 1. / CommonVars[24]*CommonVars[14]*(max(stateVars[6], 0.) - max(stateVars[16], 0.)); // A73
	dvdt[17] = 1. / CommonVars[24]*CommonVars[14]*(max(stateVars[7], 0.) - max(stateVars[17], 0.)); // A74 
	dvdt[18] = 1. / CommonVars[25]*CommonVars[15]*(max(stateVars[8], 0.) - max(stateVars[18], 0.)); // A75
	dvdt[19] = 1. / CommonVars[25]*CommonVars[15]*(max(stateVars[9], 0.) - max(stateVars[19], 0.)); // A76
	dvdt[20] = 1. / CommonVars[9]*CommonVars[16]*(max(stateVars[10], 0.) - max(stateVars[20], 0.)); // A71
	dvdt[21] = 1. / CommonVars[9]*CommonVars[16]*(max(stateVars[11], 0.) - max(stateVars[21], 0.)); // A72 
	dvdt[22] = 1. / CommonVars[10]*CommonVars[17]*(max(stateVars[12], 0.) - max(stateVars[22], 0.)); // A67
	dvdt[23] = 1. / CommonVars[10]*CommonVars[17]*(max(stateVars[13], 0.) - max(stateVars[23], 0.)); // A68
	dvdt[24] = 1. / CommonVars[11]*CommonVars[18]*(max(stateVars[14], 0.) - max(stateVars[24], 0.));// A69
	dvdt[25] = 1. / CommonVars[11]*CommonVars[18]*(max(stateVars[15], 0.) - max(stateVars[25], 0.));// A70 

	dvdt[26] = 1. / (CommonVars[12])* (CommonVars[22]*(max(stateVars[22], 0.) - max(0., stateVars[26])) + CommonVars[23]*(max(stateVars[24], 0.) - max(0., stateVars[26])) + CommonVars[21]*(max(stateVars[20], 0.) - max(0., stateVars[26])) + CommonVars[19]*(max(stateVars[16], 0.) - max(0., stateVars[26])) + CommonVars[20]*(max(stateVars[18], 0.) - max(0., stateVars[26]))); // A77
	dvdt[27] = 1. / (CommonVars[12])* (CommonVars[22]*(max(stateVars[23], 0.) - max(0., stateVars[27])) + CommonVars[23]*(max(stateVars[25], 0.) - max(0., stateVars[27])) + CommonVars[21]*(max(stateVars[21], 0.) - max(0., stateVars[27])) + CommonVars[19]*(max(stateVars[17], 0.) - max(0., stateVars[27])) + CommonVars[20]*(max(stateVars[19], 0.) - max(0., stateVars[27])));  // A78
}

void gasTransport::additionalRelations() {
	sao2 = (AlveoliLungsAlg[5] - AlveoliLungsAlg[9])/LungsGasParams[16]/1.34; // A56;
}

void gasTransport::output(){
	sampleGASstateVars << scientific << setprecision(12);
    sampleGASstateVars << time << " ";
	// State variables
	for (int i = 0; i < stateVars.size(); i++) {    
		sampleGASstateVars << stateVars[i] << " ";     
	}
	sampleGASstateVars << "\n";
	sampleGASstateVars.flush();
	// Algebraic variables
	sampleGASalveoliLungsAlg << scientific << setprecision(12);
	sampleGASalveoliLungsAlg << time << " ";
	for (int i = 0; i < AlveoliLungsAlg.size(); i++) {
		sampleGASalveoliLungsAlg << AlveoliLungsAlg[i] << " ";
	}
	sampleGASalveoliLungsAlg << pHsysArt << " ";
	sampleGASalveoliLungsAlg << pHpulCap << " ";
	sampleGASalveoliLungsAlg << "\n";
	sampleGASalveoliLungsAlg.flush();
	
	sampleGAStissueVenousGasAlg << scientific << setprecision(12);
	sampleGAStissueVenousGasAlg << time << " ";
	for (int i = 0; i < TissueVenousGasAlg.size(); i++) {
		sampleGAStissueVenousGasAlg << TissueVenousGasAlg[i] << " ";
	}
	sampleGAStissueVenousGasAlg << "\n";
	sampleGAStissueVenousGasAlg.flush();
}
