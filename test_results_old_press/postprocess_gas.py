import numpy as np
import matplotlib.pyplot as plt
import os as os

from pylab import *
from tabulate import *
    
pathToCase = "./"

def exploreTRANSPORTsimulation(path,period=None):
    
    # path output
    pathOut = path+"postproc_gas/"
    if not os.path.exists(pathOut):
            os.makedirs(pathOut)
    
    # read cvs results from files
    test = path+"GAS_"
    gasData = np.genfromtxt(test+"stateVars.txt") 
    gasDataAlg = np.genfromtxt(test+"AlveoliLungsAlg.txt") 
    
    # define final time
    tEnd = gasData[-1,0]
    print("tEnd: %.3f" % tEnd)
    
    T = period
    print("T: %.3f" % T)
    
    # get last respiratory cycle mask
    tIdx = gasData[:,0]>=(tEnd-T)
    # get time for last cycle 
    t_real = gasData[tIdx,0]
    t=np.zeros_like(t_real)
    for i in range(len(t_real)):
         t[i]=t_real[i]-t_real[0]
    time = gasDataAlg[:,0]
    # # get variables of interest
    pao2 = gasDataAlg[tIdx,8]
    paco2 = gasDataAlg[tIdx,9]
    cvo2 = gasData[tIdx,27]
    cvco2 = gasData[tIdx,28]
    pAco2 = gasDataAlg[tIdx,3]
    
    chpo2 = gasData[tIdx,13]
    chpco2 = gasData[tIdx,14]
    cbpo2 = gasData[tIdx,15]
    cbpco2 = gasData[tIdx,16]
    cspo2 = gasData[tIdx,9]
    cspco2 = gasData[tIdx,10]
    cmpo2 = gasData[tIdx,11]
    cmpco2 = gasData[tIdx,12]
    cepo2 = gasData[tIdx,7]
    cepco2 = gasData[tIdx,8]
    
    FAo2 = gasData[tIdx,3]
    FAco2 = gasData[tIdx,4]

    mdoto2 = gasDataAlg[tIdx,11]
    mdotco2 = gasDataAlg[tIdx,12]
    

    pao2avg = np.average(pao2)
    paco2avg = np.average(paco2)
    #pvo2avg = np.average(pvo2) 
    #pvco2avg = np.average(pvco2)
    cvo2avg = np.average(cvo2)
    cvco2avg = np.average(cvco2)
    #pAo2avg = np.average(pAo2)
    pAco2avg = np.average(pAco2)
    #pdo2max = np.average(pdo2)
    #pdco2min = np.average(pdco2)
    
    chpo2avg = np.average(chpo2)
    chpco2avg = np.average(chpco2)
    cbpo2avg = np.average(cbpo2)
    cbpco2avg = np.average(cbpco2)
    cspo2avg = np.average(cspo2)
    cspco2avg = np.average(cspco2)
    cmpo2avg = np.average(cmpo2)
    cmpco2avg = np.average(cmpco2)
    cepo2avg = np.average(cepo2)
    cepco2avg = np.average(cepco2)
    
    #sao2avg = np.average(sao2)
    
    FAo2avg = np.average(FAo2)
    FAco2avg = np.average(FAco2)

    f = open(pathOut+'tab10.out','w')
    f.write("pao2: %.4f mmHg\n" % pao2avg )
    f.write("paco2: %.4f #mmHg\n" % paco2avg)
    #f.write("pvo2: %.4f mmHg\n" % pvo2avg) 
    #f.write("pvco2: %.4f mmHg\n" % pvco2avg )
    f.write("cvo2: %.4f ml/ml\n" % cvo2avg )
    f.write("cvco2: %.4f ml/ml\n" % cvco2avg)
    #f.write("pAo2: %.4f mmHg\n" % pAo2avg)
    f.write("pAco2: %.4f mmHg\n" % pAco2avg )
    #f.write("pdo2: %.4f mmHg\n" % pdo2max)
    #f.write("pdco2: %.4f mmHg\n" % pdco2min)
    
    
    f.write("chpo2: %.4f ml/ml\n" % chpo2avg)
    f.write("chpco2: %.4f ml/ml\n" % chpco2avg)
    f.write("cbpo2: %.4f ml/ml\n" % cbpo2avg)
    f.write("cbpco2: %.4f ml/ml\n" % cbpco2avg)
    f.write("cspo2: %.4f ml/ml\n" % cspo2avg)
    f.write("cspco2: %.4f ml/ml\n" % cspco2avg)
    f.write("cepo2: %.4f ml/ml\n" % cepo2avg)
    f.write("cepco2: %.4f ml/ml\n" % cepco2avg)
    f.write("cmpo2: %.4f ml/ml\n" % cmpo2avg)
    f.write("cmpco2: %.4f ml/ml\n" % cmpco2avg)
    #f.write("sao2: %.4f #ml/ml\n" % sao2avg)
    f.write("FAo2: %.4f mmHg\n" % FAo2avg )
    f.write("FAco2: %.4f mmHg\n" % FAco2avg)
    
    f.close()
    
    # f = open(pathOut+'gasFAO2CO2.out','w')
    # f.write("FAo2: %.4f mmHg\n" % FAo2avg )
    # f.write("FAco2: %.4f mmHg\n" % FAco2avg)
    # f.close()
    
    # plt.figure(figsize=(10,3))
    # plt.title("Flux Mdoto2-Mdotco2")
    # plt.plot(t,mdoto2,'r-',label="$mdoto2$ [%.2f(%.2f/%.2f]" % (np.average(mdoto2),mdoto2.min(),mdoto2.max()))
    # plt.plot(t,mdotco2,'g-',label="$mdotco2$ [%.2f(%.2f/%.2f]" % (np.average(mdotco2),mdotco2.min(),mdotco2.max()))
    # plt.xlabel("$t\,[s]$")
    # plt.ylabel("$Q\,[mmHg/s]$")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(pathOut+"oxygenCarbonDioxideFlux"+figType)
    # plt.close()
    
    
    # plt.figure()
    # plt.title("Normalized state variables")
    # labelsgas = np.arange(1,gasData.shape[1]-1)
    # labelsgas = labelsgas.tolist()
    # N = time[tIdx].shape[0]
    # for i in range(gasData.shape[1]-1):
    #     av = np.convolve(gasData[:,i+1], np.ones(N)/N, mode='valid')
    #     diff = time.shape[0]-av.shape[0]
    #     plt.plot(time[diff:],av/np.average(gasData[tIdx,i+1])) 
    # plt.xlabel(r'$t\,[s]$')
    # plt.ylabel(r'$\bar{y}\,[-]$')
    # plt.tight_layout()
    # plt.legend(labelsgas)
    # plt.savefig(pathOut+"transient"+figType)
    # plt.close()
    
    plt.figure(figsize=(10,3))
    plt.title("$O_2$ partial pressure in systemic arteries, systemic veins and alveoli")
    plt.plot(t,pao2,'r-',label="$P_{a,O_2}$ [%.2f(%.2f/%.2f]" % (np.average(pao2),pao2.min(),pao2.max()))
   # plt.plot(t,pvo2,'b-',label="$P_{v,O_2}$ [%.2f(%.2f/%.2f]" % (np.average(pvo2),pvo2.min(),pvo2.max()))
   # plt.plot(t,pAo2,'g-',label="$P_{A,O_2}$ [%.2f(%.2f/%.2f]" % (np.average(pAo2),pAo2.min(),pAo2.max()))
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$P\,[mmHg]$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(pathOut+"partialpressureO2.png")
    plt.close()

    plt.figure()
    plt.title("$CO_2$ partial pressure in systemic arteries, systemic veins and alveoli")
    plt.plot(t,paco2,'r-',label="$P_{a,CO_2}$ [%.2f(%.2f/%.2f]" % (np.average(paco2),paco2.min(),paco2.max()))
    #plt.plot(t,pvco2,'b-',label="$P_{v,CO_2}$ [%.2f(%.2f/%.2f]" % (np.average(pvco2),pvco2.min(),pvco2.max()))
    plt.plot(t,pAco2,'g-',label="$P_{A,CO_2}$ [%.2f(%.2f/%.2f]" % (np.average(pAco2),pAco2.min(),pAco2.max()))
    plt.xlabel("$t\,[s]$")
    plt.ylabel("$P\,[mmHg]$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(pathOut+"partialpressureCO2.png")
    plt.close()
    
    # plt.figure(figsize=(10,3))
    # plt.title("$O_2$ and $CO_2$ concentration in mixed venous blood")
    # plt.plot(t,pvo2,'r-',label="$C_{v,O_2}$ [%.2f(%.2f/%.2f]" % (np.average(cvo2),cvo2.min(),cvo2.max()))
    # plt.plot(t,pvco2,'g-',label="$C_{v,CO_2}$ [%.2f(%.2f/%.2f]" % (np.average(cvco2),cvco2.min(),cvco2.max()))
    # plt.xlabel("$t\,[s]$")
    # plt.ylabel("$C\,[mL/mL]$")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(pathOut+"concentrationVeins.png")
    # plt.close()
    
    # plt.figure()
    # plt.title("Hemoglobin $O_2$ saturation")
    # plt.plot(t,sao2,'r-',label="$S_{a,O_2}$ [%.2f(%.2f/%.2f]" % (np.average(sao2),sao2.min(),sao2.max()))
    # plt.xlabel("$t\,[s]$")
    # plt.ylabel("$\%\,[-]$")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(pathOut+"HemoOxygenSat.png")
    # plt.close()
    



nameCVSdiff = "cvsOutput.out"
nameCVSalg = "cvsOutputAlg.out"
nameLUNGdiff = "lungOutput.out"
nameLUNGalg = "lungOutputAlg.out"
nameTRANSPORTdiff = "gasOutput.out"
nameTRANSPORTalg = "gasOutputAlg.out"
nameCONTROLdiff = "controlOutput.out"
nameCONTROLalg = "controlOutputAlg.out"
nameCONTROLRESPdiff = "controlRespOutput.out"
nameCONTROLRESPalg = "controlRespOutputAlg.out"

#savePlots(pathToCase,nameCVSdiff,nameCVSalg,period =3000.)
#exploreCVS_per_tabella(pathToCase,nameCVSdiff,nameCVSalg,period = 60.)
#exploreCVS_controlled_variables(pathToCase,nameCONTROLdiff,nameCONTROLalg,nameCVSalg,nameCVSdiff,period = 60.)
#exploreCVSsimulation(pathToCase,nameCVSdiff,nameCVSalg,period = 60.)
#exploreLUNGsimulation(pathToCase,nameLUNGdiff,nameLUNGalg,period=60.)
exploreTRANSPORTsimulation(pathToCase,period=60.)
#exploreCONTROLsimulation(pathToCase,nameCONTROLdiff,nameCONTROLalg,period=60.)
#exploreCONTROLRESPsimulation(pathToCase,nameCONTROLRESPdiff,nameCONTROLRESPalg,period=60.)





print("DONE")  
    
    
    
