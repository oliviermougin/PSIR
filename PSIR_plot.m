clear all,
% 7T tissue parameters:
% my = myelin, wm = white matter,
% csf = cerebrospinal fluid, gm = grey matter
T1my=300;T1wm=1200;T1csf=5000;T1gm=2000;
Mmy=0.71;Mwm=0.71;Mgm=0.83;Mcsf=1;


%% %% UK7T 0.7mm PSIR SSi=3.5
% Scanning parameters: for details see dx.doi.org/10.1002/mrm.26061
% SSi: shot to shot interval, TR: repetition time
% SENSE: acceleration factor
% SSi=3500;TR=6.3;SSii=10*SSi;
% SENSEy=3;SENSEz=1;
% TFE=224;
% FAt1=5;FAt2=2;TI1min=715;TI1a=725;TIbis2=5;%CNRt=1.39/NPI=-0.1
% FB=180;inv2=false;

SSi=5000;TR=16;SSii=10*SSi;
SENSEy=3;SENSEz=1;
TFE=63;
FAt1=8;FAt2=5;TI1min=570;TI1a=1080;TIbis2=150;%TIbis2=PhaseInterval-2*TI1min
FB=180;inv2=false;


% B1 inhomogeneity in percent
B1inh=100;
FA1=FAt1*B1inh/100;FA2=FAt2*B1inh/100;
TIbis = TI1a + ceil((TFE/2)*TR)+TIbis2;
% Estimated acquisition time for whole head (in min)
time=SSi*TFE/0.7/SENSEy/SENSEz/60000;

fig=true;

% Initialisation
Mmy0 = Mmy*(1-exp(-SSi/T1my));
Mwm0 = Mwm*(1-exp(-SSi/T1wm));
Mgm0 = Mgm*(1-exp(-SSi/T1gm));
Mcsf0 = Mcsf*(1-exp(-SSi/T1csf));
% Time variable
T = 1:10*SSi;
M0=ones(1,numel(T))*Mmy0;
M1=ones(1,numel(T))*Mwm0;
M2=ones(1,numel(T))*Mgm0;
M3=ones(1,numel(T))*Mcsf0;

TI1i=10*TI1a;
TIbisi=10*TIbis;
TI4i=TIbisi+TI1i; % Needed to correct for phase

% Beware, linear encoding for the PSIR:
% centre of the k-space encoded at the middle of the TFE train
% so Inversion time TI set at centre of TFE train
TFE_bis=TFE/2;
TI0i=TI1i-ceil(TFE_bis*10*TR);
TI2i=TI1i+ceil(TFE_bis*10*TR);
TI3i=TI4i-ceil(TFE_bis*10*TR);
TI5i=TI4i+ceil(TFE_bis*10*TR);

if TI0i<0
    disp('Error, wrong timing, please reajust!')
    return
end

% Necessary to repeat up to 5 times the acquisition
% to enter into a steady-state
for shot=1:5
    
    % Inversion (supposed to be perfect)
    M0(1)=cosd(FB)*M0(1);
    M1(1)=cosd(FB)*M1(1);
    M2(1)=cosd(FB)*M2(1);
    M3(1)=cosd(FB)*M3(1);
    
    % relaxation before the acquisition
    for i=2:1:TI0i
        M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
        M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
        M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
        M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
    end
    %disp(['White and grey matter abs diff before acq ' num2str(abs(M1(TI0))-abs(M2(TI0))) ])
    % start of acquisition of inv1
    for i=TI0i:1:TI1i
        if (mod(i-TI0i,10*TR)==0)
            M0(i)=Mmy+(cosd(FA1)*M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(cosd(FA1)*M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(cosd(FA1)*M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(cosd(FA1)*M3(i-1)-Mcsf)*exp(-.1/T1csf);
        else
            M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
        end
    end
    
    % TI1: acquisition of the centre of k-space
    % NPI is the contrast in the inv1 image between tissues
    CNRt = abs(M1(TI1i)-M2(TI1i)) ./ time;
    CNR = abs(M1(TI1i)-M2(TI1i));
    NPI = M1(TI1i)+M2(TI1i);
    NPI2 = M2(TI1i)+M3(TI1i);
    
    %% Signal in image 1 for WM
    INV1=M1(TI1i);
    
    for i=TI1i:1:TI2i
        if (mod(i-TI0i,10*TR)==0)
            M0(i)=Mmy+(cosd(FA1)*M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(cosd(FA1)*M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(cosd(FA1)*M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(cosd(FA1)*M3(i-1)-Mcsf)*exp(-.1/T1csf);
        else
            M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
        end
    end
    
    for i=TI2i+1:TIbisi
        M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
        M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
        M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
        M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
    end
    
    % If inversion before the second acquisition (bug in the patch)
    if inv2
        M0(TIbisi)=cosd(FB)*M0(TIbisi);
        M1(TIbisi)=cosd(FB)*M1(TIbisi);
        M2(TIbisi)=cosd(FB)*M2(TIbisi);
        M3(TIbisi)=cosd(FB)*M3(TIbisi);
    end
    for i=TIbisi+1:TI3i
        M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
        M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
        M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
        M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
    end
    
    for i=TI3i:1:TI4i
        if (mod(i-TI3i,10*TR)==0)
            M0(i)=Mmy+(cosd(FA2)*M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(cosd(FA2)*M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(cosd(FA2)*M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(cosd(FA2)*M3(i-1)-Mcsf)*exp(-.1/T1csf);
        else
            M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
        end
    end
    
    %% Signal in image 2 for WM
    INV2=M1(TI4i);
    
    for i=TI4i:1:TI5i
        if (mod(i-TI3i,10*TR)==0)
            M0(i)=Mmy+(cosd(FA2)*M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(cosd(FA2)*M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(cosd(FA2)*M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(cosd(FA2)*M3(i-1)-Mcsf)*exp(-.1/T1csf);
        else
            M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
            M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
            M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
            M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
        end
    end
    
    for i=TI5i:SSii
        M0(i)=Mmy+(M0(i-1)-Mmy)*exp(-.1/T1my);
        M1(i)=Mwm+(M1(i-1)-Mwm)*exp(-.1/T1wm);
        M2(i)=Mgm+(M2(i-1)-Mgm)*exp(-.1/T1gm);
        M3(i)=Mcsf+(M3(i-1)-Mcsf)*exp(-.1/T1csf);
    end
    
    % Reasignment of magnetization for the next shot
    M0(1)=M0(SSii);M1(1)=M1(SSii);M2(1)=M2(SSii);M3(1)=M3(SSii);
end

if fig,
    figure
    subplot(2,1,1), hold on
    plot(T./10,M0,'k')
    plot(T./10,0.5*(M1+M0),'y')
    plot(T./10,M1,'g')
    plot(T./10,0.5*(M1+M2),'m')
    plot(T./10,M2,'b')
    plot(T./10,0.5*(M2+M3),'c')
    plot(T./10,M3,'r')
    xlim([-10 SSi]);
    
    xlabel('Time (in ms)')
    ylabel('Longitudinal Magnetization (Normalised)')
    TI0 = TI0i./10*ones(1,100);
    TI1 = TI1i./10*ones(1,100);
    TI2 = TI2i./10*ones(1,100);
    TIbisn = TIbisi./10*ones(1,100);
    TI3 = TI3i./10*ones(1,100);
    TI4 = TI4i./10*ones(1,100);
    TI5 = TI5i./10*ones(1,100);
    linee = -1:2/99:1;
    line(zeros(1,100),linee,'LineStyle','-','MarkerSize',10,'Color','b')
    line(TI0,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    line(TI1,linee,'LineStyle',':','MarkerSize',3,'Color','k')
    line(TI2,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    line(TIbisn,linee,'LineStyle','-','MarkerSize',10,'Color','b')
    line(TI3,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    line(TI4,linee,'LineStyle',':','MarkerSize',3,'Color','k')
    line(TI5,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    plot(1:SSi,zeros(SSi,1),':','MarkerSize',3,'Color','k')
    legend({'Myelin water','Mix myelin+white matter','White Matter','Mix WM+GM','Grey Matter','Mix GM+CSF','CSF'},'Location','SouthEast')
    
    subplot(2,1,2), hold on
    plot(T./10,M1,'k')
    plot(T./10,abs(M1-M2),'m')
    plot(T./10,M2,'c')
    plot(T./10,abs(M2-M3),'r')
    plot(T./10,M3,'g')
    xlim([-10 SSi]);
    linee = -1:2/99:1;
    line(zeros(1,100),linee,'LineStyle','-','MarkerSize',10,'Color','b')
    line(TI0,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    line(TI1,linee,'LineStyle',':','MarkerSize',3,'Color','k')
    line(TI2,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    line(TIbisn,linee,'LineStyle','-','MarkerSize',10,'Color','b')
    line(TI3,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    line(TI4,linee,'LineStyle',':','MarkerSize',3,'Color','k')
    line(TI5,linee,'LineStyle','-','MarkerSize',3,'Color','k')
    plot(1:SSi,zeros(SSi,1),':','MarkerSize',3,'Color','k')
    legend({'WM';'WM-GM';'GM';'GM-CSF';'CSF'})
end

PSIR = INV1/(abs(INV1)+abs(INV2));
MP2RAGE = INV1*INV2/(INV1.^2+INV2.^2);

disp(['The signal in the PSIR image is ' num2str(PSIR) ', and the signal in the MP2RAGE is ' num2str(MP2RAGE)]);
disp(['The signal in the INV1 image is ' num2str(INV1) ', and the signal in the INV2 is ' num2str(INV2)]);



