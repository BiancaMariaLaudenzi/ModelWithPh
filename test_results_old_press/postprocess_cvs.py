import numpy as np
import matplotlib.pyplot as plt
import os as os
import re
from pylab import *
from scipy.signal import find_peaks
from tabulate import *

pathToCase = './' #LungsGasControls

def parse_custom_file(file_path):
    data = {}
    with open(file_path, 'r') as file:
        current_section = None
        for line in file:
            if line.strip().startswith('['):
                current_section = line.strip()[1:-1]
                data[current_section] = {}
            elif current_section is not None:
                if line.strip() and not line.strip().startswith('//'):
                    key, value = map(str.strip, line.split('='))
                    data[current_section][key] = {'value': float(value.split(';')[0])}
    return data

def get_value(param, data):
    if param == 'vusa' or  param =='vuep' or  param =='vusp' or  param =='vump' or  param =='vuhp' or  param =='vubp'or  param =='vuhv'or  param =='vubv' or param =='vupa' or  param =='vupp' or param =='vupv':
        cathegory='unstressed_volumes'
    if param =='tvvmax'or  param == 'tvvu' or  param == 'tvvmin':
        cathegory='thoracicVeins'
    if param=='vulv' or  param ==  'vurv' or param=='p0lv' or param=='kelv' or param== 'kerv':
        cathegory='ventricles'
    if param =='vula' or  param == 'vura' or param=='tcLA':
        cathegory='atria'
    if param == 'vsa'  or  param == 'vep' or  param == 'vsp' or  param == 'vmp' or  param == 'vhp' or  param == 'vbp' or  param == 'vev' or  param == 'vsvnew' or  param == 'vmvnew' or  param == 'vhv' or  param == 'vbv' or  param == 'vtv'  or  param == 'vpa' or  param ==  'vpp' or  param ==  'vpv' or  param == 'vps' or  param == 'vla' or  param == 'vlv'  or  param == 'vra'  or  param == 'vrv':
        cathegory='initial_conditions'
    if param=='k_closeMV':
        cathegory='valve'
    if param=='T0':
        cathegory='basalHR'
    if param=='GTv':
        cathegory='heartControlParams'
    if param=='pn':
        cathegory='baroAfferentParam'
    if param=='vusv0' or param== 'vuev0' or param=='vumv0':
        cathegory='cardiovascularControl'   
    if param == 'csatco2' or param =='h2' or param == 'k2':
        cathegory='LungsGas'
    if param in data[cathegory]:
        value=data[cathegory][param]['value']
    return value

def exploreCVS_per_tabella(pathDir,period = 60):
    
    # path output
    pathOut = pathDir+"postproc_CVS/"
    if not os.path.exists(pathOut):
            os.makedirs(pathOut)

    path=pathDir+'CVS_'
    cvsHeart = np.genfromtxt(path+'Heart.out')
    cvsPulCirc = np.genfromtxt(path+'PulCirc.out')
    cvsStateVars = np.genfromtxt(path+'stateVars.out')
    cvsSysCirc = np.genfromtxt(path+'SysCirc.out')
    cvsValve = np.genfromtxt(path+'Valve.out')
    path_params='../baseline/'
    cvsParams = parse_custom_file(path_params+'/cardiovascular.dat')
    commonParams = parse_custom_file(path_params+'/commonParams.dat')
    # define final time
    tEnd = cvsHeart[-1,0]
    print("tEnd: %.3f" % tEnd)
    
    # define last cardiac cycle period
    T = period
    print("T: %.3f" % T)
   
    # get last cardiac cycle mask
    tIdx = (cvsHeart[:,0]>=(tEnd-T))&(cvsHeart[:,0]<=(tEnd)) #vettore di veri e falsi
    t = cvsHeart[tIdx,0]
    
    u = cvsHeart[tIdx,2]
    

    nSample = 200
    [end,ini] = getPeaks(t,u,nSample)
    t_ini=t[ini]
    # print(t)
    #print(t_ini)
    #print(len(t_ini))
    tIdx = cvsHeart[:,0]>=t_ini[0] #vettore di veri e falsi
    t = cvsHeart[tIdx,0]
    vlv = cvsStateVars[tIdx,18]
    psa = cvsSysCirc[tIdx,1]
    [vlv_max,vlv_min] = getPeaks(t,vlv,nSample)
    [psa_max,psa_min] = getPeaks(t,psa,nSample)
    LVEDV = vlv[vlv_max]
    LVESV = vlv[vlv_min]
    SBP = psa[psa_max]
    DBP = psa[psa_min]
    LVESP= psa[vlv_min]
    BSA=1.92 #1.6488 #m^2
    psaavg=np.zeros(len(t_ini)-1)
    praavg=np.zeros(len(t_ini)-1)
    ptvavg=np.zeros(len(t_ini)-1)
    ppvavg=np.zeros(len(t_ini)-1)
    ppaavg=np.zeros(len(t_ini)-1)
    qAVavg=np.zeros(len(t_ini)-1)
    qtvavg=np.zeros(len(t_ini)-1)
    qbpavg=np.zeros(len(t_ini)-1)
    qhpavg=np.zeros(len(t_ini)-1)
    PSBP=np.zeros(len(t_ini)-1)
    PDBP=np.zeros(len(t_ini)-1)
    PP_A=np.zeros(len(t_ini)-1)
    PPP_A=np.zeros(len(t_ini)-1)
    LVSV=np.zeros(len(t_ini)-1)
    LVEF=np.zeros(len(t_ini)-1)
    ESBP1= np.zeros(len(t_ini)-1)
    ESBP2= np.zeros(len(t_ini)-1)
    E_LV=np.zeros(len(t_ini)-1)
    E_a=np.zeros(len(t_ini)-1)
    ratio=np.zeros(len(t_ini)-1)
    E_LV2=np.zeros(len(t_ini)-1)
    E_a2=np.zeros(len(t_ini)-1)
    ratio2=np.zeros(len(t_ini)-1)
    CO=np.zeros(len(t_ini)-1)
    SV_fromCO=np.zeros(len(t_ini)-1)
    Csa=np.zeros(len(t_ini)-1)
    qAVavg_I=np.zeros(len(t_ini)-1)
    CO_I=np.zeros(len(t_ini)-1)
    vtotavg=np.zeros(len(t_ini)-1)
    vtot_unstressed_avg=np.zeros(len(t_ini)-1)
    vHeartAvg=np.zeros(len(t_ini)-1)
    vPulAvg=np.zeros(len(t_ini)-1)
    vVenAvg=np.zeros(len(t_ini)-1)
    vArtAvg=np.zeros(len(t_ini)-1)
    psvavg = np.zeros(len(t_ini)-1)
    pevavg = np.zeros(len(t_ini)-1)
    pmvavg = np.zeros(len(t_ini)-1)
    phvavg = np.zeros(len(t_ini)-1)
    pbvavg = np.zeros(len(t_ini)-1)

    for i in range(len(t_ini)-1):
        t_diff=t_ini[i+1]-t_ini[i]
        #print(t_ini[i])
        # get each cardiac cycle mask
        a=cvsHeart[:,0]>=t_ini[i]
        b=cvsHeart[:,0]<=t_ini[i+1] #vettore di veri e falsi
        # print(a)
        # print(b)
        tIdx=(a & b)
        #print(tIdx)
        # get time for each cycle 
        t = cvsHeart[tIdx,0]
        psa = cvsSysCirc[tIdx,1]
        ptv = cvsSysCirc[tIdx,14]
        ppv = cvsPulCirc[tIdx,2]
        ppa = cvsPulCirc[tIdx,1]
        pra = cvsHeart[tIdx,8]
        qAV = cvsStateVars[tIdx,28]
        qtv = cvsSysCirc[tIdx,26]
        qbp = cvsSysCirc[tIdx,12]
        qhp = cvsSysCirc[tIdx,11]
        vlv = cvsStateVars[tIdx,18]
        vla = cvsStateVars[tIdx,17]
        vrv = cvsStateVars[tIdx,21]
        vra = cvsStateVars[tIdx,20]
        vsa = cvsStateVars[tIdx,1]
        vsp = cvsStateVars[tIdx,3]
        vep = cvsStateVars[tIdx,2]
        vmp = cvsStateVars[tIdx,4]
        vhp = cvsStateVars[tIdx,5]
        vbp = cvsStateVars[tIdx,6]
        vev = cvsStateVars[tIdx,7]
        vsv = cvsStateVars[tIdx,8]
        vmv = cvsStateVars[tIdx,9]
        vhv = cvsStateVars[tIdx,10]
        vbv = cvsStateVars[tIdx,11]
        vtv = cvsStateVars[tIdx,12]
        vpa = cvsStateVars[tIdx,13]
        vpp = cvsStateVars[tIdx,14]
        vps = cvsStateVars[tIdx,16]
        vpv = cvsStateVars[tIdx,15]

        vusa = get_value('vusa', cvsParams)
        vuep = get_value('vuep', cvsParams)
        vump = get_value('vump', cvsParams)
        vuhp = get_value('vuhp', cvsParams)
        vubp = get_value('vubp', cvsParams)
        vusp = get_value('vusp', cvsParams)
        vuhv = get_value('vuhv', cvsParams)
        vubv = get_value('vubv', cvsParams)
        tvvu = get_value('tvvu', cvsParams)
        vupa = get_value('vupa', cvsParams)
        vupp = get_value('vupp', cvsParams)
        vupv = get_value('vupv', cvsParams)
        vula = get_value('vula', cvsParams)
        vura = get_value('vura', cvsParams)
        vulv = get_value('vulv', cvsParams)
        vurv = get_value('vurv', cvsParams)
        if pathToCase == "./doCardioLungsGasControls/":
            path_control = 'CardioControl_'
            controlAlg = np.genfromtxt(pathToCase+path_control+'algebraicVars.txt')
            vumv = controlAlg[tIdx,38]
            vusv = controlAlg[tIdx,39]
            vuev = controlAlg[tIdx,40]
        else:
            vusv = get_value('vusv0', commonParams)
            vuev = get_value('vuev0', commonParams)
            vumv = get_value('vumv0', commonParams)

        vtot = vsa + \
                vsp + vep + vmp + vhp + vbp +\
                vev + vhv + vbv + vmv + vsv +\
                vtv +\
                vpa + vpp + vps + vpv +\
                vra + vrv + vla + vlv
        vtot_unstressed = vusa +\
                        vusp + vuep + vump + vuhp + vubp +\
                        vuhv + vubv + vuev + vumv+ vusv +\
                        tvvu +\
                        vupa + vupp + vupv +\
                        vula + vura  + vulv + vurv        
        V_pul = vpa + vpp + vps + vpv 
        V_heart =vra + vrv + vla + vlv
        V_ven =vev + vhv + vbv + vmv + vsv + vtv
        V_art =vsa + vsp + vep + vmp + vhp + vbp
        psv = cvsSysCirc[tIdx,17]
        pev = cvsSysCirc[tIdx,16]
        pmv = cvsSysCirc[tIdx,18]
        phv = cvsSysCirc[tIdx,19]
        pbv = cvsSysCirc[tIdx,20]
        
        psaavg[i] = np.average(psa)    #media su ogni ciclo cardiaco
        ptvavg[i] = np.average(ptv)
        ppvavg[i] = np.average(ppv)
        ppaavg[i] = np.average(ppa)
        praavg[i] = np.average(pra) 
        qAVavg[i] = np.average(qAV)
        qAVavg_I[i] = np.average(qAV/1000*60/BSA)
        qtvavg[i] = np.average(qtv)
        qbpavg[i] = np.average(qbp)
        qhpavg[i] = np.average(qhp)
        PSBP[i] = 120/108*SBP[i]
        PDBP[i] = 74/75*DBP[i]
        PP_A[i]= SBP[i]-DBP[i]
        PPP_A[i]= PSBP[i]-PDBP[i]
        LVSV[i]=LVEDV[i]-LVESV[i]
        LVEF[i]=LVSV[i]/LVEDV[i]
        ESBP1[i] = (2*PSBP[i]+PDBP[i])/3
        ESBP2[i] = 120/108*LVESP[i]
        E_LV[i]=ESBP1[i]/LVESV[i]*BSA
        E_a[i]=ESBP1[i]/LVSV[i]*BSA
        ratio[i]=E_a[i]/E_LV[i]
        E_LV2[i]=ESBP2[i]/LVESV[i]*BSA
        E_a2[i]=ESBP2[i]/LVSV[i]*BSA
        ratio2[i]=E_a2[i]/E_LV2[i]
        CO[i]=LVSV[i]/t_diff
        SV_fromCO[i]=np.average(qAV)*t_diff
        CO_I[i]=LVSV[i]/t_diff/1000*60/BSA
        Csa[i]=LVSV[i]/PPP_A[i]/BSA
        vtotavg[i] = np.average(vtot)
        vtot_unstressed_avg[i] = np.average(vtot_unstressed)
        vHeartAvg[i] = np.average(V_heart)
        vPulAvg[i] = np.average(V_pul)
        vVenAvg[i] = np.average(V_ven)
        vArtAvg[i] = np.average(V_art)
        psvavg[i] = np.average(psv)
        pevavg[i] = np.average(pev)
        pmvavg[i] = np.average(pmv)
        phvavg[i] = np.average(phv)
        pbvavg[i] = np.average(pbv)
            
    qAVstd = np.std(qAVavg)  #std sul vettore delle medie di ogni ciclo cardiaco
    qAVstd_I = np.std(qAVavg_I) 
    qtvstd = np.std(qtvavg)
    psastd = np.std(psaavg)
    ptvstd = np.std(ptvavg)
    ppvstd = np.std(ppvavg)
    prastd = np.std(praavg)
    qbpstd = np.std(qbpavg)
    qhpstd = np.std(qhpavg)
    ppastd = np.std(ppaavg) 
    LVEDVstd = np.std(LVEDV)
    LVESVstd = np.std(LVESV)
    DBPstd = np.std(DBP)
    SBPstd= np.std(SBP)
    PDBPstd = np.std(PDBP)
    PSBPstd= np.std(PSBP)
    PP_Astd = np.std(PP_A) 
    PPP_Astd = np.std(PPP_A)
    LVSVstd = np.std(LVSV)
    SV_fromCOstd = np.std(SV_fromCO)
    LVEFstd=np.std(LVEF)
    E_LVstd=np.std(E_LV)
    E_a_std=np.std(E_a)
    ratioStd=np.std(ratio)
    E_LVstd2=np.std(E_LV2)
    E_a_std2=np.std(E_a2)
    ratioStd2=np.std(ratio2)
    COstd=np.std(CO)
    COstd_I=np.std(CO_I)
    CsaStd=np.std(Csa)
    vtotStd=np.std(vtot)
    vtot_unstressedStd=np.std(vtot_unstressed)
    vArtStd=np.std(V_art)
    vPulStd=np.std(V_pul)
    vHeartStd=np.std(V_heart)
    vVenStd=np.std(V_ven)
    psvstd = np.std(psvavg)
    pevstd = np.std(pevavg)
    pmvstd = np.std(pmvavg)
    phvstd = np.std(phvavg)
    pbvstd = np.std(pbvavg)
        

    psaavg = np.average(psaavg)    #media delle medie su ogni ciclo cardiaco
    ptvavg = np.average(ptvavg)
    ppvavg = np.average(ppvavg)
    ppaavg = np.average(ppaavg)
    praavg = np.average(praavg)
    qAVavg = np.average(qAVavg)
    qAVavg_I = np.average(qAVavg_I)
    qtvavg = np.average(qtvavg)
    qbpavg = np.average(qbpavg)
    qhpavg = np.average(qhpavg)
    LVEDVavg = np.average(LVEDV)
    LVESVavg = np.average(LVESV)
    DBPavg = np.average(DBP)
    SBPavg = np.average(SBP)
    PDBPavg = np.average(PDBP)
    PSBPavg = np.average(PSBP)
    PP_Aavg = np.average(PP_A)
    PPP_Aavg = np.average(PPP_A)
    LVSVavg = np.average(LVSV)
    SV_fromCOavg = np.average(SV_fromCO)
    LVEFavg=np.average(LVEF)
    E_LVavg=np.average(E_LV)
    E_a_avg=np.average(E_a)   
    ratioAvg=np.average(ratio)
    E_LVavg2=np.average(E_LV2)
    E_a_avg2=np.average(E_a2)   
    ratioAvg2=np.average(ratio2)
    COavg=np.average(CO)
    COavg_I=np.average(CO_I)
    CsaAvg=np.average(Csa)
    vtotavg=np.average(vtotavg)
    vtot_unstressedAvg=np.average(vtot_unstressed_avg)
    vArtAvg=np.average(V_art)
    vPulAvg=np.average(V_pul)
    vHeartAvg=np.average(V_heart)
    vVenAvg=np.average(V_ven)
    psvavg = np.average(psvavg)
    pevavg = np.average(pevavg)
    pmvavg = np.average(pmvavg)
    phvavg = np.average(phvavg)
    pbvavg = np.average(pbvavg)
    HR = len(t_ini)/T*60


    f=open(pathOut+'Table_Haemodynamics_variables_T='+str(period)+'s_BSA='+str(BSA)+'.out','w')
    table=[['\hline'],['\\textbf{Parameter}','\\textbf{Units}','\\textbf{Model}','\\textbf{Ref. val.}','\\textbf{ Ref.}'],['\hline'],['\hline'],\
           ['$\mathrm{HR}$', 'beats/min', '%.1f'% (HR),  '68(12)' ,'\cite{mceniery2008central}' ],\
            ['$\mathrm{LVSV}$', 'mL', '%.4f [%.4f]' % (LVSVavg,LVSVstd), '112(19) ',  '\cite{hudsmith2005normal}'],\
            ['$\mathrm{LVSV}$ from $\mathrm{CO}$', 'mL', '%.4f [%.4f]' % (SV_fromCOavg,SV_fromCOstd), '112(19) ',  '\cite{hudsmith2005normal}'],\
            ['$\mathrm{LVSV}$ index', 'mL/m$^2$', '%.4f [%.4f]' % (LVSVavg/BSA, LVSVstd), '56(8) ', '\cite{hudsmith2005normal}'],\
            ['$\mathrm{V_{LV,max}}$', 'mL', '%.4f [%.4f]' % (LVEDVavg, LVEDVstd),'160(29) ','\cite{hudsmith2005normal}'],\
            ['$\mathrm{LVEF}$',  '%'  , '%.4f [%.4f]' % (LVEFavg, LVEFstd), '63(5)' ,'\cite{borlaug2009contractility}'],\
            ['$C_{sa}\\mathrm{I}$', 'mL/mmHg/m$^2$' , '%.4f [%.4f]' % (CsaAvg, CsaStd), '1.08(0.27)', '\cite{abdelhammed2005noninvasive}'],\
            ['\hline'],\
            ['\\textcolor{blue}{Formula from \cite{najjarAgeGenderAffect2004}}'],\
            ['\hline'],\
            ['$\mathrm{E}_\mathrm{LV}\mathrm{I}$', 'mmHg$\cdot$m$^2$/mL',  '%.4f [%.4f]' % (E_LVavg, E_LVstd) , '4.50(-)'  ,'\cite{najjarAgeGenderAffect2004}'],\
            ['$\mathrm{Ea}\mathrm{I}$' , 'mmHg$\cdot$m$^2$/mL',  '%.4f [%.4f]' % (E_a_avg, E_a_std) , '2.20(-)' ,'\cite{najjarAgeGenderAffect2004}'],\
            ['$\mathrm{Ea}\mathrm{I}/\mathrm{E}_\mathrm{LV}\mathrm{I}$', '-' , '%.4f [%.4f]' % (ratioAvg, ratioStd), '0.58(-)' , '\cite{najjarAgeGenderAffect2004}'],\
            ['\hline'],\
            ['$\\textcolor{red}{[P_{sa}[V_{LV,min}]]_{Perif}}^{*}$ '],\
            ['\hline'],\
            ['$\mathrm{E}_\mathrm{LV}\mathrm{I}$', 'mmHg$\cdot$m$^2$/mL', '%.4f [%.4f]' % (E_LVavg2, E_LVstd2), '4.50(-)' ,'\cite{najjarAgeGenderAffect2004}'],\
            ['$\mathrm{Ea}\mathrm{I}$','mmHg$\cdot$m$^2$/mL',  '%.4f [%.4f]' % (E_a_avg2, E_a_std2) ,'2.20(-)' ,'\cite{najjarAgeGenderAffect2004}'],\
            ['$\mathrm{Ea}\mathrm{I}/\mathrm{E}_\mathrm{LV}\mathrm{I}$' ,'-' , '%.4f [%.4f]' % (ratioAvg2, ratioStd2)  , '0.58(-)', '\cite{najjarAgeGenderAffect2004}'],\
            ['\hline'],['\hline'],\
            ['\\textbf{Pressure}','\\textbf{Units}','\\textbf{Model}', '\\textbf{Ref. val.}','\\textbf{ Ref.}' ],\
            ['\hline'],['\hline'],\
            ['$\mathrm{MAP}$'   , 'mmHg','%.4f [%.4f]' % (psaavg, psastd), '88(8)' ,'\cite{mceniery2008central} '],\
            ['$\mathrm{CDBP}$' ,  'mmHg','%.4f [%.4f]' % (DBPavg, DBPstd), '74(9)','\cite{mceniery2008central} '],\
            ['$\mathrm{CSBP}$' , 'mmHg' ,'%.4f [%.4f]' % (SBPavg, SBPstd), '103(8)' ,'\cite{mceniery2008central} '],\
            ['$\mathrm{PDBP}$' , 'mmHg' ,'%.4f [%.4f]$^{*}$' % (PDBPavg, PDBPstd) ,  '73(8)','\cite{mceniery2008central} '],\
            ['$\mathrm{PSBP}$' , 'mmHg' ,'%.4f [%.4f]$^{*}$' % (PSBPavg, PSBPstd) ,'123(10)','\cite{mceniery2008central} '],\
            ['$\mathrm{MPAP}$'  ,  'mmHg' ,'%.4f [%.4f]' % (ppaavg, ppastd) , '14(3)','\cite{lauRestingPulmonaryArtery}'],\
            ['$\mathrm{CVP}$ ' , 'mmHg' ,'%.4f [%.4f]' % (ptvavg, ptvstd),'(0-5)' ,' \cite{levickIntroductionCardiovascularPhysiology2010}' ],\
            ['\hline'],\
            ['$\mathrm{CPP}$' , 'mmHg' ,'%.4f [%.4f]' % (PP_Aavg, PP_Astd), '29(5)' ,' \cite{mceniery2008central} ' ],\
            ['$\mathrm{PPP}$', 'mmHg','%.4f [%.4f]' % (PPP_Aavg, PPP_Astd), '50(9)' ,'\cite{mceniery2008central} ' ],\
            ['\hline'],['\hline'],\
            ['\\textbf{Cardiac cycle average flow}','\\textbf{Units}','\\textbf{Model}', '\\textbf{Ref. val.}','\\textbf{ Ref.}' ],\
            ['\hline'],['\hline'],\
            ['Cardiac output (qAV)', 'mL/s','%.4f [%.4f]' % (qAVavg, qAVstd) , '83.3(33.3)' ,'\cite{cattermoleNormalRangesCardiovascular2017}'],\
            ['Cardiac output (LVSV$\\times$t$_{\\text{cycle}}$)', 'mL/s' ,'%.4f [%.4f]' % (COavg, COstd) ,'83.3(33.3)' ,'\cite{cattermoleNormalRangesCardiovascular2017}'],\
            ['CO index (qAV)', 'mL/s/m$^2$','%.4f [%.4f]' % (qAVavg_I, qAVstd_I),'2.9(0.8)','\cite{ganau1992patterns}'],\
            ['CO index (LVSV$\\times$t$_{\\text{cycle}}$)', 'mL/s/m$^2$','%.4f [%.4f]' %(COavg_I, COstd_I),'2.9(0.8)','\cite{ganau1992patterns}'],\
            ['Venous flow ', 'mL/s','%.4f [%.4f]' % (qtvavg, qtvstd) , '---' ,'---'],\
            ['Cerebral blood flow' ,'mL/s','%.4f [%.4f]' %  (qbpavg, qbpstd) ,'12.18(2.12)' ,'\cite{fordCharacterizationVolumetricFlow2005a}'],\
            ['Coronary blood flow' ,'mL/s ,'  ,'%.4f [%.4f]' % (qhpavg, qhpstd),'4.5(1.36)' ,'\cite{sakamotoRelationDistributionCoronary2013}'],['\hline'],['\hline']]

    t=tabulate(table, tablefmt='latex')
    s=convert2LatexString(t)
    f.write(s)
    # save result
    f=open(pathOut+'Haemodynamics_variables_T='+str(period)+'s.out','w')
    f.write("V_tot: %.4f [%.4f] #mL\n" % (vtotavg, vtotStd) )
    f.write("V_tot_unstressed: %.4f [%.4f] #mL\n" % (vtot_unstressedAvg, vtot_unstressedStd) )
    f.write("HR: %.4f #beats/min\n" % HR )
    f.write("LVSV: %.4f [%.4f] #mL\n" % (LVSVavg, LVSVstd) )
    f.write("LVEF: %.4f [%.4f] \n" % (LVEFavg, LVEFstd) )
    f.write("C_sa: %.4f [%.4f] \n" % (CsaAvg, CsaStd) )
    f.write("LVEDV: %.4f [%.4f] #mL\n" % (LVEDVavg, LVEDVstd) )
    f.write("LVESV: %.4f [%.4f] #mL\n" % (LVESVavg, LVESVstd) )  
    f.write("LVSV_I: %.4f [%.4f] #mL\n" % (LVSVavg/BSA, LVSVstd) )
    f.write("#------------------------------------------\n")
    f.write("#con formula \n")
    f.write("E_LV: %.4f [%.4f] \n" % (E_LVavg/BSA, E_LVstd) ) 
    f.write("E_a: %.4f [%.4f] \n" % (E_a_avg/BSA, E_a_std) )
    f.write("Ratio: %.4f [%.4f] \n" % (ratioAvg, ratioStd) )
    f.write("E_LV_I: %.4f [%.4f] \n" % (E_LVavg, E_LVstd) ) 
    f.write("E_a_I: %.4f [%.4f] \n" % (E_a_avg, E_a_std) )
    f.write("#------------------------------------------\n")
    f.write("#con Psa[min_Vlv] \n") 
    f.write("E_LV: %.4f [%.4f] \n" % (E_LVavg2/BSA, E_LVstd2) ) 
    f.write("E_a: %.4f [%.4f] \n" % (E_a_avg2/BSA, E_a_std2) ) 
    f.write("Ratio: %.4f [%.4f] \n" % (ratioAvg2, ratioStd2) )
    f.write("E_LV_I: %.4f [%.4f] \n" % (E_LVavg2, E_LVstd2) ) 
    f.write("E_a_I: %.4f [%.4f] \n" % (E_a_avg2, E_a_std2) )

    f.write("#------------------------------------------\n")
    f.write("#Pressures \n")
    f.write("#------------------------------------------\n")
    f.write("psa: %.4f [%.4f] #mmHg\n" % (psaavg, psastd) )
    f.write("CDBP: %.4f [%.4f] #mmHg\n" % (DBPavg, DBPstd) )
    f.write("CSBP: %.4f [%.4f] #mmHg\n" % (SBPavg, SBPstd) ) 
    f.write("PDBP: %.4f [%.4f] #mmHg\n" % (PDBPavg, PDBPstd) )
    f.write("PSBP: %.4f [%.4f] #mmHg\n" % (PSBPavg, PSBPstd) )
    f.write("ppa: %.4f [%.4f] #mmHg\n" % (ppaavg, ppastd) )
    f.write("ptv: %.4f [%.4f] #mmHg\n" % (ptvavg, ptvstd) )
    f.write("ppv: %.4f [%.4f] #mmHg\n" % (ppvavg, ppvstd) )
    f.write("RAP: %.4f [%.4f] #mmHg\n" % (praavg, prastd) )
    f.write("#------------------------------------------\n")
    f.write("CPP_A: %.4f [%.4f] #mmHg\n" % (PP_Aavg, PP_Astd) ) 
    f.write("PPP_A: %.4f [%.4f] #mmHg\n" % (PPP_Aavg, PPP_Astd) )
    f.write("#------------------------------------------\n")
    f.write("psv: %.4f [%.4f] #mL\n" % (psvavg, psvstd) )
    f.write("pev: %.4f [%.4f] #mL\n" % (pevavg, pevstd) )
    f.write("pmv: %.4f [%.4f] #mL\n" % (pmvavg, pmvstd) )
    f.write("phv: %.4f [%.4f] #mL\n" % (phvavg, phvstd) )
    f.write("pbv: %.4f [%.4f] #mL\n" % (pbvavg, pbvstd) )

    f.write("#------------------------------------------\n")
    f.write("#Blood flow distribution \n")
    f.write("#------------------------------------------\n")
    f.write("qAV: %.4f [%.4f] #mL/s\n" % (qAVavg, qAVstd) )
    f.write("CO: %.4f [%.4f] #mL/s \n" % (COavg, COstd) )
    f.write("qtv: %.4f [%.4f] #mL/s\n" % (qtvavg, qtvstd) )
    f.write("qbp: %.4f [%.4f] #mL/s\n" % (qbpavg, qbpstd) )
    f.write("qhp: %.4f [%.4f] #mL/s\n" % (qhpavg, qhpstd) )
    f.write("CI_qAV: %.4f [%.4f] #L/min/m^2\n" % (qAVavg_I, qAVstd_I) )
    f.write("CI_CO: %.4f [%.4f] #L/min/m^2\n" % (COavg_I, COstd_I) )

    plt.close()
    print('Heamo vars: DONE')

    def func_simul(pct, allvalues):
        absolute = int(pct / 100.*np.sum(allvalues))
        return f"{np.round(pct,2)}%\n({absolute} mL)"
    figure, axis = plt.subplots(nrows=1, ncols=2, figsize = (10,5))
    figure.suptitle(f'Total Volume ({int(vtotavg)} mL)', fontweight='bold')
    mylabels = ['Stressed volume [20 - 25 %]', "Unstressed volume [75 - 80 %]"] #from Davis
    volumes = [vtotavg-vtot_unstressedAvg, vtot_unstressedAvg] 
    mylabels2 = ['Pulmonary Circulation [9 %]', "Heart [7 %]", 'Venous Blood [69 %]', 'Arterial Blood [15 %]'] # from Celant-Muller - I added capillary blood volume [5%] in venous volume [64%]
    volumes2 = [vPulAvg , vHeartAvg, vVenAvg, vArtAvg] 
    pie = axis[0].pie(volumes, labels = mylabels, autopct = lambda pct: func_simul(pct, volumes), shadow=True, startangle=90)
    pie2 = axis[1].pie(volumes2, labels = mylabels2, autopct = lambda pct: func_simul(pct, volumes2), shadow=True, startangle=90)
    for label in pie[1]:
        label.set_fontsize(8)
    for label in pie2[1]:
        label.set_fontsize(8)
    figure.tight_layout()
    figure.savefig(pathOut+"volume_pie.pdf")


def convert2LatexString(s):
    r = r'(\^\{\})'; s = re.sub(r, "^", s)
    s = re.sub(r'\\([\$\_\{\}\^])', r'\1', s)
    s = re.sub(r'(\\textbackslash{})', r'\\', s)
    return s

def getPeaks(t,d,nSample):
    done = 0
    height = 0.9
    distance = 50
#while done==0:
    plt.plot(t,d)
    peaks_max, _ = find_peaks(d, height=height*np.max(d),distance=distance)
    peaks_min, _ = find_peaks(-d, height=-height*np.max(d),distance=distance)
    plt.plot(t[peaks_max], d[peaks_max], "x")
    plt.plot(t[peaks_min], d[peaks_min], "o")
    #plt.show() #block=False

    # val = input("Are you satisfied with peak selection [y/n]: ")
    # if val=='y':
    #     done = 1
    # else :
    #     height = float(input("Insert new height [%.2f]: " % height))
    #     distance = int(input("Insert new distance [%d]: " % distance))
    #     print(height,distance)  
    plt.close()
    return [peaks_max,peaks_min]   

def savePlots(path,pathDiff,pathAlg,figType='.png',period = 60.):

    pathOut = path+"postproc/Plots/"
    if not os.path.exists(pathOut):
            os.makedirs(pathOut)
    
    # define template cvs model to access variables
    cvs = cvsClass.cardioVascularSystem("templateForPostprocessing/inputFiles/cardiovascular.dat",outFolder='templateForPostprocessing/')
    lung = lungMechClass.lungMechanics("templateForPostprocessing/inputFiles/lungMechanics.dat",outFolder='templateForPostprocessing/')
    control = cardiovascularcontrolClass.cardioVascularControl("templateForPostprocessing/inputFiles/control.dat",outFolder='templateForPostprocessing/')
    controlResp = respiratoryControlClass.respiratoryControl("templateForPostprocessing/inputFiles/controlResp.dat",outFolder='templateForPostprocessing/')
    
    # read cvs results from files
    cvsData = np.genfromtxt(path+pathDiff) # these are cvs model state variables
    cvsDataAlg = np.genfromtxt(path+pathAlg) # these are "algebraic" variables
    # read cvs results from files
    lungDataAlg = np.genfromtxt(path+"lungOutputAlg.out")
    controlDataAlg = np.genfromtxt(path+"controlOutputAlg.out")
    controlData = np.genfromtxt(path+"controlOutput.out")
    controlRespDataAlg = np.genfromtxt(path+"controlRespOutputAlg.out")
    #print(cvsData)
    n1 = cvsData.shape[0]
    n2 = cvsDataAlg.shape[0]
    n = min(n1,n2)
    #print(n)
    #print(cvsData)
    cvsData = cvsData[:n-20,:]
    cvsDataAlg = cvsDataAlg[:n-20,:]
    controlData = controlData[:n-20,:]
    controlDataAlg = controlDataAlg[:n-20,:]
    lungDataAlg = lungDataAlg[:n-20,:]
    controlRespDataAlg = controlRespDataAlg[:n-20,:]
    # define final time
    tEnd = cvsData[-1,0]
    print("tEnd: %.3f" % tEnd)
    
    # define last cardiac cycle period
    if period==None:
            # get last cardiac cycle mask
            T = 2.5
            tIdx = cvsDataAlg[:,0]>=(tEnd-T)
            t = cvsDataAlg[tIdx,0]
            u = cvsDataAlg[tIdx,cvs.xAlgNames.index('u')+1]
            nSample = 200
            [end,ini] = getPeaks(t,u,nSample)
            t_ini=t[ini]
            #t_end=t[end]
            #print(t_ini)
            a=cvsDataAlg[:,0]>=t_ini[0]
            b=cvsDataAlg[:,0]<=t_ini[1] #vettore di veri e falsi
            # print(a)
            # print(b)
            tIdx=(a & b)
    else:
            T = period
            # get last cardiac cycle mask
            tIdx = cvsDataAlg[:,0]>=(tEnd-T)
            
    print("T: %.3f" % T)
    t_real = cvsDataAlg[tIdx,0]
    t=np.zeros_like(t_real)
    for i in range(len(t)):
         t[i]=t_real[i]-t_real[0]
    rr=controlRespDataAlg[tIdx,controlResp.xAlgNames.index('rr')+1]
    HP=controlDataAlg[tIdx,control.xAlgNames.index('T')+1]
    ptau=controlData[tIdx,control.xNames.index('ptau')+1]
    fab=controlDataAlg[tIdx,control.xAlgNames.index('fab')+1]
    #print(self.cvs.xAlgNames.index('qAV')+1)
    # xiAV=cvsData[t0:,self.cvs.xNames.index('xiAV')+1]
    # xiPV=cvsData[t0:,self.cvs.xNames.index('xiPV')+1]
    qAV=cvsData[tIdx,cvs.xNames.index('qAV')+1]
    qMV=cvsData[tIdx,cvs.xNames.index('qMV')+1]
    
    ppl=lungDataAlg[tIdx,lung.xAlgNames.index('ppl1')+1] #pressure pleura
    plv=cvsDataAlg[tIdx,cvs.xAlgNames.index('plv')+1]
    psa=cvsDataAlg[tIdx,cvs.xAlgNames.index('psa')+1]
    pla=cvsDataAlg[tIdx,cvs.xAlgNames.index('pla')+1]
    pmaxlv=cvsDataAlg[tIdx,cvs.xAlgNames.index('pmaxlv')+1]
    emaxlv=cvsDataAlg[tIdx,cvs.xAlgNames.index('emaxlv')+1]
    #xiPV=cvsData[t0,self.cvs.xNames.index('xiPV')+1]
    #xiPV=cvsData[t0,self.cvs.xNames.index('xiPV')+1]
    qPV=cvsData[tIdx,cvs.xNames.index('qPV')+1]
    qTV=cvsData[tIdx,cvs.xNames.index('qTV')+1]
    
    prv=cvsDataAlg[tIdx,cvs.xAlgNames.index('prv')+1]
    ppa=cvsDataAlg[tIdx,cvs.xAlgNames.index('ppa')+1]
    pra=cvsDataAlg[tIdx,cvs.xAlgNames.index('pra')+1]
    pmaxrv=cvsDataAlg[tIdx,cvs.xAlgNames.index('pmaxrv')+1]
    
    phi=cvsDataAlg[tIdx,cvs.xAlgNames.index('phi')+1]
    vlv=cvsData[tIdx,cvs.xNames.index('vlv')+1]
    vrv=cvsData[tIdx,cvs.xNames.index('vrv')+1]
      
    vla = cvsData[tIdx,cvs.xNames.index('vla')+1]
    vra = cvsData[tIdx,cvs.xNames.index('vra')+1]

    vtv = cvsData[tIdx,cvs.xNames.index('vtv')+1]
    ptvtm = cvsDataAlg[tIdx,cvs.xAlgNames.index('ptvtm')+1]
    
    eLA=cvsDataAlg[tIdx,cvs.xAlgNames.index('eLA')+1]
    #print(eLA)
    eRA=cvsDataAlg[tIdx,cvs.xAlgNames.index('eRA')+1]
    phi=cvsDataAlg[tIdx,cvs.xAlgNames.index('phi')+1]

    vpa = cvsData[tIdx,cvs.xNames.index('vpa')+1]
    vpp = cvsData[tIdx,cvs.xNames.index('vpp')+1]
    vpv = cvsData[tIdx,cvs.xNames.index('vpv')+1]

    sys=phi*emaxlv*(vlv-14.7576)
    dia=(1.-phi)*1.5*(np.exp(0.014*vlv)-1.)
    sys_hyp=phi*1.3*emaxlv*(vlv-14.7576)
    dia_hyp=(1.-phi)*1.3*1.5*(np.exp(0.014*vlv)-1.)
    dia_new=(1.-phi)*1.5*(np.exp(0.014*1.3*vlv)-1.)

      
    # fig = plt.figure(dpi=400)
    
    # ax1 = fig.add_subplot(321)
    # l1, = ax1.plot(t, qAV, 'm--',label='qAV') 
    # v1, = ax1.plot(t, qMV, 'y--',label='qMV') 
    # ax1.grid()
    # ax1.set_ylabel('q [mL/s]')
    # ax1.set_xticklabels([])
      
    # ax5 = fig.add_subplot(322)
    # l2,=ax5.plot(t, qPV, 'm-',label='qPV') 
    # v2, = ax5.plot(t, qTV, 'y-',label='qTV') 
    # ax5.grid()
    # ax5.set_xticklabels([])

    
    # ax2 = fig.add_subplot(323)
    # l3,=ax2.plot(t, plv, 'r--',label='plv') 
    # l4,=ax2.plot(t, pla, 'b--',label='pla') 
    # l5,=ax2.plot(t, psa, 'g--',label='psa') 
    # ax2.set_ylabel('p [mmHg]')
    # ax2.grid() 
    # ax2.set_xticklabels([])
    
    # ax4 = fig.add_subplot(324)
    # l6,=ax4.plot(t, prv, 'r-',label='prv') 
    # l7,=ax4.plot(t, pra, 'b-',label='pra') 
    # l8,=ax4.plot(t, ppa, 'g-',label='ppa') 
    # ax4.grid() 
    # ax4.set_xticklabels([])
  
      
    # ax6 = fig.add_subplot(325)
    # l9,=ax6.plot(t, vlv, 'c--',label='vlv') 
    # l10,=ax6.plot(t, vla, 'k--',label='vla')
    # ax6.set_ylabel('V [mL]')
    # ax6.set_xlabel('t [s]')
    # ax6.grid()
  
    # ax8 = fig.add_subplot(326)
    # l11,=ax8.plot(t, vrv, 'c-',label='vrv') 
    # l12,=ax8.plot(t, vra, 'k-',label='vra') 
    # ax8.set_xlabel('t [s]')
    # ax8.grid()

    # leg1=plt.legend([l1,v1,l3,l4,l5,l9,l10],['qAV [%.2f(%.2f/%.2f)]' % (np.average(qAV),qAV.min(),qAV.max()),'qMV [%.2f(%.2f/%.2f)]' % (np.average(qMV),qMV.min(),qMV.max()),'plv [%.2f(%.2f/%.2f)]' % (np.average(plv),plv.min(),plv.max()),'pla [%.2f(%.2f/%.2f)]' % (np.average(pla),pla.min(),pla.max()),'psa [%.2f(%.2f/%.2f)]' % (np.average(psa),psa.min(),psa.max()),'vlv [%.2f(%.2f/%.2f)]' % (np.average(vlv),vlv.min(),vlv.max()),'vla [%.2f(%.2f/%.2f)]' % (np.average(vla),vla.min(),vla.max())],bbox_to_anchor=(-0.1,-0.5),title = "Left Heart")
    # leg2=plt.legend([l2,v2,l6,l7,l8,l11,l12],['qPV [%.2f(%.2f/%.2f)]' % (np.average(qPV),qPV.min(),qPV.max()),'qTV [%.2f(%.2f/%.2f)]' % (np.average(qTV),qTV.min(),qTV.max()),'prv [%.2f(%.2f/%.2f)]' % (np.average(prv),prv.min(),prv.max()),'pra [%.2f(%.2f/%.2f)]' % (np.average(pra),pra.min(),pra.max()),'ppa [%.2f(%.2f/%.2f)]' % (np.average(ppa),ppa.min(),ppa.max()),'vrv [%.2f(%.2f/%.2f)]' % (np.average(vrv),vrv.min(),vrv.max()),'vra [%.2f(%.2f/%.2f)]' % (np.average(vra),vra.min(),vra.max())],bbox_to_anchor=(1.1,-0.5),title = "Right Heart")
    # gca().add_artist(leg1)
    # #plt.legend()
    # plt.tight_layout() 
    # plt.savefig(pathOut+"Heart.png")
    # #plt.show()
    # plt.close() 
    
    # fig = plt.figure(dpi=400)
    # ax1 = fig.add_subplot(111)
    # ax1.plot(t, eLA, 'r-',label='eLA (Atrii)') 
    # ax1.plot(t, eRA, 'g-',label='eRA (Atrii)') 
    # ax1.plot(t, phi, 'b-',label='phi (Ventricle)')
    # ax1.legend()
    # plt.savefig(pathOut+"Elastances.png")
    # #plt.show()
    # plt.close()

    # # left ventricle P-V loop
    
    # def poly_area2D(poly):
    #     """
        
    #     """
    #     total = 0.0
    #     N = len(poly)
    #     for i in range(N):
    #         v1 = poly[i]
    #         v2 = poly[(i+1) % N]
    #         total += v1[0]*v2[1] - v1[1]*v2[0]
    #     return abs(total/2)
    
    # p =np.zeros((vlv.shape[0]+1,2))
    # p[:-1,0] = vlv
    # p[-1,0] = vlv[0]
    # p[:-1,1] = plv
    # p[-1,1] = plv[0]
    # #print(p)
    # area = poly_area2D(p)
    # work = area*133.3322*1e-6
    # plt.title("P-V curve for left ventricle")
    # plt.plot(vlv,plv,'r-',label="$W_{lv}$ = %.3f J" % (work))
    # plt.ylabel("$P\,[mmHg]$")
    # plt.xlabel("$V\,[mL]$")
    # plt.legend()
    # plt.savefig(pathOut+"pvLeftVentricle"+figType)
    # plt.close()
    # # pressure in thorax # Eq. 2
    # tvkxp = 2. # mmHg
    # tvvu = 130. # ml
    # tvkxv = 8. # ml
    # tvd1 = 0.3855 # mmHg
    # tvk1 = 0.15 # mmHg/ml
    # tvkr = 0.001 # mmHg s / ml
    # tvvmax = 350. # ml
    # tvvu = 130. # ml
    # tvr0 = 0.025 # mmHg s / ml
    # tvd2 = -5. # mmHg
    # tvk2 = 0.4 # mmHg
    # tvvmin = 50. # ml
    # tvkxp = 2. # mmHg
    # tvkxv = 8. # ml
    # vtv_tot=np.linspace(-50,100)
    # psi=np.zeros_like(vtv_tot)
    # ptvtm_tot=np.zeros_like(vtv_tot)
    # for i in range(len(vtv_tot)):
    #     psi[i] =  tvkxp/np.exp((( vtv_tot[i]+ tvvu)/ tvkxv)-1.)
    #     if ( vtv_tot[i]+ tvvu)>= tvvu:
    #         ptvtm_tot[i] =  tvd1 +  tvk1*(( tvvu+ vtv_tot[i])- tvvu) - psi[i]
    #     else:
    #         ptvtm_tot[i] =  tvd2 +  tvk2*np.exp(( tvvu+ vtv_tot[i])/ tvvmin)-psi[i]

    # plt.title("P-V curve for thoracic vein")
    # plt.plot(vtv_tot,ptvtm_tot,'b-')
    # plt.plot(vtv,ptvtm,'r-',label='P [%.2f(%.2f/%.2f]' % (np.average(ptvtm),ptvtm.min(),ptvtm.max()))
    # plt.plot(vtv,ptvtm,'r-',label='V [%.2f(%.2f/%.2f]'% (np.average(vtv),vtv.min(),vtv.max()))
    # plt.ylabel("$P\,[mmHg]$")
    # plt.xlabel("$V\,[mL]$")
    # plt.legend()
    # plt.savefig(pathOut+"pvThoracicVein"+figType)
    # plt.close()

    # plt.title("Pressure LV")
    # plt.plot(t,ppl,'r-',label='P$_{pl}$ [%.2f(%.2f/%.2f)](%.2f)' % (np.average(ppl),ppl.min(),ppl.max(),ppl[0]))
    # plt.plot(t,plv,'b-',label='P$_{LV}$ [%.2f(%.2f/%.2f)](%.2f)'% (np.average(plv),plv.min(),plv.max(),plv[0]))
    # plt.ylabel("$P\,[mmHg]$")
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.savefig(pathOut+"Pressure LV"+figType)
    # plt.close()


    # fig = plt.figure(dpi=400)
    # ax1 = fig.add_subplot(111)
    # #ax1.plot(t, pmaxlv, 'r-',label='pLVmax [%.2f(%.2f/%.2f]' % (np.average(pmaxlv),pmaxlv.min(),pmaxlv.max()))
    # ax1.plot(t, sys, 'g--',label='sys [%.2f(%.2f/%.2f]' % (np.average(sys),sys.min(),sys.max()))
    # ax1.plot(t, dia, 'b--',label='dia [%.2f(%.2f/%.2f]' % (np.average(dia),dia.min(),dia.max()))
    # ax1.plot(t, sys_hyp, 'r--',label='sys_hyp [%.2f(%.2f/%.2f]' % (np.average(sys_hyp),sys_hyp.min(),sys_hyp.max()))
    # ax1.plot(t, dia_hyp, 'y--',label='dia_hyp [%.2f(%.2f/%.2f]' % (np.average(dia_hyp),dia_hyp.min(),dia_hyp.max()))
    # #ax1.plot(t, dia_new, 'y--',label='dia_new [%.2f(%.2f/%.2f]' % (np.average(dia_new),dia_new.min(),dia_new.max()))
    # ax1.legend()
    # plt.savefig(pathOut+"pLVmax.png")
    # #plt.show()
    # plt.close()
    
    # elv=plv/vlv
    # elvmax=pmaxlv/vlv
    # nSample = 200
    # [max_e,min_e] = getPeaks(t,elvmax,nSample)
    # t_min=t[min_e]
    # a= cvsDataAlg[:,0]>=t_min[0] #vettore di veri e falsi
    # b=cvsDataAlg[:,0]<=t_ini[1] #vettore di veri e falsi
    # tIdx_e =(a & b)
    # t_e = cvsDataAlg[tIdx_e,0]
    # pmaxlv_b = cvsDataAlg[tIdx_e,cvs.xAlgNames.index('pmaxlv')+1]
    # vlv_b = cvsData[tIdx_e,cvs.xNames.index('vlv')+1]
    # e_b=pmaxlv_b/vlv_b
    # fig = plt.figure(dpi=400)
    # ax1 = fig.add_subplot(111)
    # #ax1.plot(t, pmaxlv, 'r-',label='pLVmax [%.2f(%.2f/%.2f]' % (np.average(pmaxlv),pmaxlv.min(),pmaxlv.max()))
    # ax1.plot(t, elv, 'r--',label='elv [%.2f(%.2f/%.2f)]' % (np.average(elv),elv.min(),elv.max()))
    # ax1.plot(t, elvmax, 'g--',label='elv_max [%.2f(%.2f/%.2f)], E_B [%.2f]' % (np.average(elvmax),elvmax.min(),elvmax.max(),np.average(e_b)))
    # ax1.legend()
    # plt.savefig(pathOut+"Elastances_new.png")
    # plt.show()
    # plt.close()
    
    # fig = plt.figure(dpi=400)
    # ax1 = fig.add_subplot(111)
    # ax1.plot(t_e, e_b, 'r--',label='E_b [%.2f(%.2f/%.2f)]' % (np.average(e_b),e_b.min(),e_b.max()))
    # ax1.legend()
    # plt.savefig(pathOut+"Passive_Elastance.png")
    # #plt.show()
    # plt.close()

    # PPGh=cvsData[tIdx,cvs.xNames.index('vhp')+1]
    # PPGb=cvsData[tIdx,cvs.xNames.index('vbp')+1]
    # PPGm=cvsData[tIdx,cvs.xNames.index('vmp')+1]
    # PPGs=cvsData[tIdx,cvs.xNames.index('vsp')+1]
    # PPGe=cvsData[tIdx,cvs.xNames.index('vep')+1]
    # PPGsa=cvsData[tIdx,cvs.xNames.index('vsa')+1]
    # Ph=cvsDataAlg[tIdx,cvs.xAlgNames.index('php')+1]
    # Pb=cvsDataAlg[tIdx,cvs.xAlgNames.index('pbp')+1]
    # Pm=cvsDataAlg[tIdx,cvs.xAlgNames.index('pmp')+1]
    # Ps=cvsDataAlg[tIdx,cvs.xAlgNames.index('psp')+1]
    # Pe=cvsDataAlg[tIdx,cvs.xAlgNames.index('pep')+1]
    # Psa=cvsDataAlg[tIdx,cvs.xAlgNames.index('psa')+1]
    # qsa=cvsDataAlg[tIdx,cvs.xAlgNames.index('qsa')+1]

    # plt.title("Peripheral PPG")
    # plt.plot(t,PPGh,'r-',label='PPG$_{h}$ [%.2f(%.2f/%.2f)]' % (np.average(PPGh),PPGh.min(),PPGh.max()))
    # plt.plot(t,PPGb,'b-',label='PPG$_{b}$ [%.2f(%.2f/%.2f)]' % (np.average(PPGb),PPGb.min(),PPGb.max()))
    # plt.plot(t,PPGm,'g-',label='PPG$_{m}$ [%.2f(%.2f/%.2f)]' % (np.average(PPGm),PPGm.min(),PPGm.max()))
    # plt.plot(t,PPGs,'c-',label='PPG$_{s}$ [%.2f(%.2f/%.2f)]' % (np.average(PPGs),PPGs.min(),PPGs.max()))
    # plt.plot(t,PPGe,'m-',label='PPG$_{e}$ [%.2f(%.2f/%.2f)]' % (np.average(PPGe),PPGe.min(),PPGe.max()))
    # plt.plot(t,PPGsa,'y-',label='PPG$_{sa}$ [%.2f(%.2f/%.2f)]' % (np.average(PPGsa),PPGsa.min(),PPGsa.max()))
    # plt.ylabel("$PPG\,[mL]$")
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.savefig(pathOut+"Peripheral PPG"+figType)
    # plt.close()

    # plt.title("Peripheral pressures")
    # plt.plot(t,Ph,'r-',label='P$_{h}$ [%.2f(%.2f/%.2f)]' % (np.average(Ph),Ph.min(),Ph.max()))
    # plt.plot(t,Pb,'b-',label='P$_{b}$ [%.2f(%.2f/%.2f)]' % (np.average(Pb),Pb.min(),Pb.max()))
    # plt.plot(t,Pm,'g-',label='P$_{m}$ [%.2f(%.2f/%.2f)]' % (np.average(Pm),Pm.min(),Pm.max()))
    # plt.plot(t,Ps,'c-',label='P$_{s}$ [%.2f(%.2f/%.2f)]' % (np.average(Ps),Ps.min(),Ps.max()))
    # plt.plot(t,Pe,'m-',label='P$_{e}$ [%.2f(%.2f/%.2f)]' % (np.average(Pe),Pe.min(),Pe.max()))
    # plt.plot(t,Psa,'y-',label='P$_{sa}$ [%.2f(%.2f/%.2f)]' % (np.average(Psa),Psa.min(),Psa.max()))
    # plt.ylabel("$P\,[mmHg]$")
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.savefig(pathOut+"Peripheral pressure"+figType)
    # plt.close()

    # plt.title("Peripheral PPG normalized")
    # plt.plot(t,(PPGh-PPGh.min())/(PPGh.max()-PPGh.min()),'r-',label='PPG$_{h}$')
    # plt.plot(t,(PPGb-PPGb.min())/(PPGb.max()-PPGb.min()),'b-',label='PPG$_{b}$')
    # plt.plot(t,(PPGm-PPGm.min())/(PPGm.max()-PPGm.min()),'g-',label='PPG$_{m}$')
    # plt.plot(t,(PPGs-PPGs.min())/(PPGs.max()-PPGs.min()),'c-',label='PPG$_{s}$')
    # plt.plot(t,(PPGe-PPGe.min())/(PPGe.max()-PPGe.min()),'m-',label='PPG$_{e}$')
    # plt.ylabel("$PPG\,[norm.]$")
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.savefig(pathOut+"Peripheral PPG normalized"+figType)
    # plt.close()

    # plt.title("Peripheral PPG - pressure normalized")
    # plt.plot(t,(Ph-Ph.min())/(Ph.max()-Ph.min()),label='P$_{h}$')
    # plt.plot(t,(Pb-Pb.min())/(Pb.max()-Pb.min()),label='P$_{b}$')
    # plt.plot(t,(Pm-Pm.min())/(Pm.max()-Pm.min()),label='P$_{m}$')
    # plt.plot(t,(Ps-Ps.min())/(Ps.max()-Ps.min()),label='P$_{s}$')
    # plt.plot(t,(Pe-Pe.min())/(Pe.max()-Pe.min()),label='P$_{e}$')
    # plt.plot(t,(Psa-Psa.min())/(Psa.max()-Psa.min()),label='P$_{sa}$')
    # #plt.plot(t,(qsa-qsa.min())/(qsa.max()-qsa.min()),label='Q$_{sa}$')
    # #plt.plot(t,(plv-plv.min())/(plv.max()-plv.min()),label='P$_{lv}$')

    # plt.plot(t,(PPGh-PPGh.min())/(PPGh.max()-PPGh.min()),label='PPG$_{h}$')
    # plt.plot(t,(PPGb-PPGb.min())/(PPGb.max()-PPGb.min()),label='PPG$_{b}$')
    # plt.plot(t,(PPGm-PPGm.min())/(PPGm.max()-PPGm.min()),label='PPG$_{m}$')
    # plt.plot(t,(PPGs-PPGs.min())/(PPGs.max()-PPGs.min()),label='PPG$_{s}$')
    # plt.plot(t,(PPGe-PPGe.min())/(PPGe.max()-PPGe.min()),label='PPG$_{e}$')
    # plt.plot(t,(PPGsa-PPGsa.min())/(PPGsa.max()-PPGsa.min()),label='PPG$_{sa}$')
    # plt.ylabel("Pulse waves [norm.]")
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.savefig(pathOut+"Peripheral PPG - pressure normalized"+figType)
    # plt.close()

    
    # plt.title("Aortic and mitral valve flow")
    # plt.plot(t, qAV, 'r-',label='$q_{AV}$') 
    # plt.plot(t, qMV, 'b-',label='$q_{MV}$') 
    # plt.ylabel('q [mL/s]')
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(pathOut+"LeftHeart_q"+figType)
    # plt.close()

      
    # plt.title("Left atrium, ventricle and systemic pressures")
    # plt.plot(t, plv, 'r-',label='$P_{LV}$') 
    # plt.plot(t, pla, 'b-',label='$P_{LA}$') 
    # plt.plot(t, psa, 'g-',label='$P_{sa}$') 
    # plt.ylabel('P [mmHg]')
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(pathOut+"LeftHeart_p"+figType)
    # plt.close()


    # plt.title("Pulmonary arteries, pheripheral and veins volumes")
    # plt.plot(t, vpa, 'r-',label='$V_{pa}$') 
    # plt.plot(t, vpp, 'g-',label='$V_{pp}$') 
    # plt.plot(t, vpv, 'b-',label='$V_{pv}$') 
    # plt.ylabel('V [mL]')
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(pathOut+"Pulmonary_circulation"+figType)
    # plt.close()
  
    # plt.title("Baroreceptors firing rate")
    # plt.plot(t, fab, 'r-',label='$f_{ab}$ [%.2f(%.2f/%.2f)]' % (np.average(fab),fab.min(),fab.max()))
    # plt.ylabel('$f_{ab}\,[spikes/s]$')
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.savefig(pathOut+"cardiovascular_control.png")
    # plt.close()

    plt.title("Respiratory rate")
    plt.plot(t, rr, 'r-',label='$RR$ [%.2f(%.2f/%.2f)]' % (np.average(rr),rr.min(),rr.max()))
    plt.ylabel('$RR\,[breaths/min]$')
    plt.xlabel("$t\,[s]$")
    plt.legend()
    plt.savefig(pathOut+"respiratory_rate.png")
    plt.close()

    # plt.title("Cardiac cycle")
    # plt.plot(t, HP, 'r-',label='$T$ [%.2f(%.2f/%.2f)]' % (np.average(HP),HP.min(),HP.max()))
    # plt.ylabel('$T\,[s]$')
    # plt.xlabel("$t\,[s]$")
    # plt.legend()
    # plt.savefig(pathOut+"period.png")
    # plt.close()
    print("Plots: DONE")





# nameCVSdiff = "cvsOutput.out"
# nameCVSalg = "cvsOutputAlg.out"
# nameLUNGdiff = "lungOutput.out"
# nameLUNGalg = "lungOutputAlg.out"
# nameTRANSPORTdiff = "gasOutput.out"
# nameTRANSPORTalg = "gasOutputAlg.out"
# nameCONTROLdiff = "controlOutput.out"
# nameCONTROLalg = "controlOutputAlg.out"
# nameCONTROLRESPdiff = "controlRespOutput.out"
# nameCONTROLRESPalg = "controlRespOutputAlg.out"

#savePlots(pathToCase,nameCVSdiff,nameCVSalg,period =3000.)
exploreCVS_per_tabella(pathToCase,period = 60.)
# exploreCVS_controlled_variables(pathToCase,nameCONTROLdiff,nameCONTROLalg,nameCVSalg,nameCVSdiff,period = 60.)
# exploreCVSsimulation(pathToCase,nameCVSdiff,nameCVSalg,period = 60.)
# exploreLUNGsimulation(pathToCase,nameLUNGdiff,nameLUNGalg,period=60.)
# exploreTRANSPORTsimulation(pathToCase,nameTRANSPORTdiff,nameTRANSPORTalg,period=60.)
# exploreCONTROLsimulation(pathToCase,nameCONTROLdiff,nameCONTROLalg,period=60.)
# exploreCONTROLRESPsimulation(pathToCase,nameCONTROLRESPdiff,nameCONTROLRESPalg,period=60.)
    
