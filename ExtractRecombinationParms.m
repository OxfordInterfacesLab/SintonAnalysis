function SurfParmCell=ExtractRecombinationParms(DeltaN,TauEff,W,Resistivity, Sitype, PlotsActive,DeltaNkey, SampleName,SurfParmCell,SRHfitOptions)
% This Function extracts all the recombination parameters from measurements
% of effective lifetime using a Sinton Instrument, and plots these
% parameters as a funciont of minority carrier density

% Input variables
% DeltaN: Column vector with the carrier excess density [cm-3]
% TauEff: Column vector with effective lifetime data corresponding to
% excess density in DeltaN in [s]
% W: wafer thickness in [cm]
% Resistivity in ohm.cm
% Si type: either 'n' or 'p'
% Plots Active: Set to 1 for plotting of variables
% DeltaNkey: Minority carrier density at which Seff and iVoc are calculated [cm-3]
% 
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultAxesFontSize',16);    
set(groot,'defaultLineLineWidth',1);
set(groot,'defaultAxesLineWidth',1);


q=1.6e-19;
Vt=(1.38e-23*300)/q; 

load('SiliconData2')
if strcmp(Sitype,'n'); ResTable=nres; nieffTable=nieffN;diffusionTable=diffusionN;
elseif strcmp(Sitype,'p');ResTable=pres; nieffTable=nieffP;diffusionTable=diffusionP;
else error('Type of substrate not well defined');end

Ndop = interp1(ResTable(:,2),ResTable(:,1),Resistivity,'linear','extrap');

%extract edges of measurement as tend to have artifacts
ReIndex=length(DeltaN);
DeltaN=DeltaN(floor(ReIndex*0.02):ceil(ReIndex*0.9));
TauEff=TauEff(floor(ReIndex*0.02):ceil(ReIndex*0.9));

%interpolate to filter for noise and accelerate fits
DeltaN2=logspace(log10(min(DeltaN)),log10(max(DeltaN)),60)';
TauEff2=interp1(DeltaN,TauEff,DeltaN2,'pchip','extrap');
DeltaN3=logspace(log10(min(DeltaN2)),log10(max(DeltaN2)),30)';
TauEff2=interp1(DeltaN2,TauEff2,DeltaN3,'pchip','extrap');
DeltaN2=DeltaN3;

%calculate the intrinsic lifetime of this speciment using Richter's param
[TauRad,TauAuger]=IntBulkLifetime(Ndop,DeltaN2);

%find the eff intrinsic carrier density ni as function of injection
dummy=permute(nieffTable,[1 3 2]);
[X,Y]=meshgrid(nieffTable(2:13,1,1),dummy(1,:,1));
[Xq,Yq]=meshgrid(DeltaN2,Resistivity);
nieff2=interp2(X,Y,dummy(2:13,:,2)',Xq,Yq,'linear')';

%find resistiivyt and amb diff coefficient as function of injection
dummy=permute(diffusionTable,[1 3 2]);
[X,Y]=meshgrid(diffusionTable(2:202,1,1),dummy(1,:,1));
[Xq,Yq]=meshgrid(DeltaN2,Resistivity);
Damb2=interp2(X,Y,dummy(2:202,:,7)',Xq,Yq,'linear')';

% SRV with infinite SRH lifetime
Seff=sqrt(Damb2.*(1./TauEff2-1./TauRad-1./TauAuger)).*tan((W/2)*sqrt((1./Damb2).*(1./TauEff2-1./TauRad-1./TauAuger)));

% Approx 6 Mackel pip.2167
Joe6=W*q*nieff2.^2.*gradient(1./TauEff2-1./TauRad-1./TauAuger)./gradient(DeltaN2);

% least squares method
options = optimoptions(@lsqnonlin,'Display','off','MaxFunctionEvaluations',1e6,...
    'FunctionTolerance',1e-12,'MaxIterations',1e3,'StepTolerance',1e-12,...
    'FiniteDifferenceType','central');
% for the SRH calculation 

%SRH fitting variables TauN, TauP, Et-Ev
SRH_0=max(1./(1./TauEff2-1./TauRad-1./TauAuger)).*[1 1 0]+[1 1 0.56];%starting point
SRH_lb=(max(1./(1./TauEff2-1./TauRad-1./TauAuger)))*[0 0 0];%lower boundary
SRH_ub=[1,1,1.05];

if SRHfitOptions(1)~=0;SRH_0(1)= SRHfitOptions(1);
SRH_lb(1)= SRHfitOptions(1);SRH_ub(1)=SRHfitOptions(1);end

if SRHfitOptions(2)~=0;SRH_0(2)= SRHfitOptions(2);
SRH_lb(2)= SRHfitOptions(2);SRH_ub(2)=SRHfitOptions(2);end

if SRHfitOptions(3)~=0; SRH_0(3)= SRHfitOptions(1);
SRH_lb(3)= SRHfitOptions(3);SRH_ub(3)=SRHfitOptions(3);end



switch Sitype
    case 'p'; p_0=Ndop;n_0=nieff2.^2./p_0;%1/cm3 equilibrium electrons concentration 
    case 'n'; n_0=-Ndop;p_0=nieff2.^2./n_0;%1/cm3 equilibrium electrons concentration 
end

errorCode=0;
% Joe calculated from Kimmerle SOLMAT 142 (2015) 116–122
for i=1:4% number of times to iterate the algorithm
    
    %calculate J0 with intrinsic bulk lifetime
    Joe2=q*gradient(nieff2.^2.*sqrt(Damb2.*(1./TauEff2-1./TauRad-1./TauAuger-1./(SRHfunction(SRH_0(1),SRH_0(2),SRH_0(3),DeltaN2,nieff2,n_0,p_0)))).*...
            tan((W/2)*sqrt((1./Damb2).*(1./TauEff2-1./TauRad-1./TauAuger-1./(SRHfunction(SRH_0(1),SRH_0(2),SRH_0(3),DeltaN2,nieff2,n_0,p_0))))))./gradient(DeltaN2);

    Joe2_pos=Joe2(Joe2>0);% only account for positive values of Joe
    DeltaN_pos=DeltaN2(Joe2>0);
    if SRHfitOptions(4)==0
        
        if isempty(DeltaN_pos)% when all J0 calculated is negative the method won't work. revert to simple 
            errorCode=1;
            errordlg('Surface recombinaiton is so high J_0 must be calculated from Seff',' Error');
            Joe2=sqrt(Damb2.*(1./TauEff2-1./TauRad-1./TauAuger-1./SRH_0)).*tan((W/2)*sqrt((1./Damb2).*(1./TauEff2-1./TauRad-1./TauAuger-1./SRH_0))).*(q*nieff2.^2)./((abs(Ndop)+DeltaN2));
            Joe2_pos=Joe2(Joe2>0);% only account for positive values of Joe
            DeltaN_pos=DeltaN2(Joe2>0);
        end

        %only analise deta for deltaN>0.5Ndop
        [~,IndexMinDeltaN]=min(abs(DeltaN_pos-0.5*abs(Ndop)));

        %calculate the gradient of J0 to find flatests section
        GradJoe2=abs(gradient(Joe2_pos(IndexMinDeltaN:end))./gradient(DeltaN_pos(IndexMinDeltaN:end)));
        [~,indexJoe]=min(GradJoe2);
        indexJoe=indexJoe+IndexMinDeltaN;

        if indexJoe<0.1*length(DeltaN_pos); indexJoe=ceil(0.1*length(DeltaN_pos)) ;
        elseif indexJoe>0.9*length(DeltaN_pos); indexJoe=ceil(0.9*length(DeltaN_pos) );
        end
        %take the flates section and make it J0 calculated
        JoeAverage=mean((Joe2_pos(ceil(indexJoe*0.9):floor(indexJoe*1.1))));

    else
        indexJoe=1;
%         Joe2=ones(length(DeltaN2),1)*SRHfitOptions(4);
        JoeAverage=SRHfitOptions(4);
    end

    

    %x=[TauSRH_n,TauSRH_p, Et-Ev]
    ObjectiveFun2 = @(x)(JoeAverage*(abs(Ndop)+DeltaN2)./(q*nieff2.^2)-...
        (sqrt(Damb2.*(1./TauEff2-1./TauRad-1./TauAuger-1./(SRHfunction(x(1),x(2),x(3),DeltaN2,nieff2,n_0,p_0)))).*...
        tan((W/2)*sqrt((1./Damb2).*(1./TauEff2-1./TauRad-1./TauAuger-1./(SRHfunction(x(1),x(2),x(3),DeltaN2,nieff2,n_0,p_0)))))));
    
    [SRH_I,resnorm]= lsqnonlin(ObjectiveFun2, SRH_0 ,SRH_lb,SRH_ub,options);
    
%     dummy1=fsolve(ObjectiveFun2,SRH_0);
    SRH_0=SRH_I;
    
    if errorCode==1; break; end %when S is so high these methods is inaccurate
    
end

SRH_0=SRHfunction(SRH_I(1),SRH_I(2),SRH_I(3),DeltaN2,nieff2,n_0,p_0);

iVoc=Vt*log((DeltaN2.*(abs(Ndop)+DeltaN2))./nieff2.^2);

[~,indexDeltaN]=min(abs(DeltaN2-DeltaNkey));


% SRV with estimated SRH lifetime using algorithm
Seff2=sqrt(Damb2.*(1./TauEff2-1./TauRad-1./TauAuger-1./SRH_0)).*tan((W/2)*sqrt((1./Damb2).*(1./TauEff2-1./TauRad-1./TauAuger-1./SRH_0)));
Seff3=JoeAverage*(abs(Ndop)+DeltaN2)./(q*nieff2.^2);
Jeo1=Seff.*(q*nieff2.^2)./((abs(Ndop)+DeltaN2));

%%

if PlotsActive
    figure('OuterPosition', [50 50 1300 1000],'Name', SampleName);
    
    subplot(2,2,1);
    loglog(DeltaN,TauEff,'o',DeltaN2,1./(1./TauRad+1./TauAuger),'-', ...
        DeltaN2, 1./(1./TauRad+1./TauAuger+(2*JoeAverage*(abs(Ndop)+DeltaN2)./(W*q*nieff2.^2))+1./SRH_0),'-',...%DeltaN,W./(2*Seff2),'-',...
        DeltaN2,SRH_0,'-',...
        DeltaN2,1./(2*JoeAverage*(abs(Ndop)+DeltaN2)./(W*q*nieff2.^2)),'-',...
        [abs(Ndop) abs(Ndop)],[min(TauEff) max(TauEff)],'--' );
    
    legend({'Exp Data','\tau_{intrinsic}','\tau_{best fit}','\tau_{SRH}','\tau_{SurfaceJ_{0S}}','N_{doping}'},'location','best');
    xlabel('\Deltan [cm^{-3}]');ylabel('\tau (s)');
    xlim([1e13 1e17]);
   
    subplot(2,2,2);
    loglog(DeltaN2,Seff,'o',DeltaN2,Seff2,'-', DeltaN2,Seff3,'.-');
    legend({'S_{eff0} from \tau_{eff} (\tau_{SRH}=\infty)',sprintf('S_{eff1} when finite SRH'),'S_{eff2} from J_{0s-avg}'},'location','best');
    xlabel('\Deltan (cm^{-3})');ylabel('S_{eff} (cm/s)');
    
    subplot(2,2,3);
    loglog(DeltaN2,1e15*(Joe6),'o',DeltaN2,1e15*Joe2,'^',...
        DeltaN2,1e15*Jeo1,[1e13 1e17],1e15*[JoeAverage JoeAverage]);
    legend({'J_{0s-6} Mackel','Kimmerle','J_{0s} from S_{eff1}','J_{0s} Avg'},'location','best');
    xlabel('\Deltan (cm^{-3})');ylabel('J_{0s} (fA/cm^2)');
    
    subplot(2,2,4);semilogx(DeltaN2,iVoc,'o')
    xlabel('\Deltan (cm^{-3})');ylabel('iV_{oc} (V)');  
end

% sprintf('Function Results:\nNdop:%1.2e \nTauEff:%1.2e\nJoeAverage:%1.2e\nDeltaN_{Joe}:%1.2e\niVoc:%1.2e\nSeff:%1.2e or %1.2e\nTauSRH:%1.2e',...
%     Ndop,TauEff(indexDeltaN),JoeAverage,DeltaN_pos(indexJoe,1), iVoc(indexDeltaN), Seff(indexDeltaN),Seff2(indexDeltaN),TauSRH)

%% Asign value to logging cell

DataIndex=size(SurfParmCell);
SurfParmCell(DataIndex(1)+1,1:DataIndex(2))={SampleName,[DeltaN2,Seff,Seff2,Seff3],[DeltaN2,Joe6,Joe2],Ndop,TauEff(indexDeltaN),JoeAverage,DeltaN_pos(indexJoe,1), iVoc(indexDeltaN), Seff(indexDeltaN),Seff2(indexDeltaN),SRH_I,resnorm};

end

function F=SRHfunction(taun,taup,Et_Ev,DeltaN,ni,n_0,p_0)
q=1.6e-19;
Vt=(1.38e-23*300)/q; 

    Et=-(1.1/2-Et_Ev);
    n_1=ni.*exp(Et/Vt);
    p_1=ni.*exp(-Et/Vt);
    F=((taup*(n_0+n_1+DeltaN)+taun*(p_0+p_1+DeltaN))./(p_0+n_0+DeltaN));

end
function [t_Aug,t_Rad]=IntBulkLifetime(Ndop,Delta_n)

q=1.60e-19;%C electron charge
k=1.38e-23;%J/K boltazmann constant
T=300;%K Temperature
ni=9.7e9;%1/cm3 Semicon intrinsic carrier consentration
if Ndop>0;p_0=Ndop;n_0=ni^2/p_0;% 
else n_0=-Ndop; p_0=ni^2/n_0;
end


%excess concentration in stable region
n_d=n_0+Delta_n;
p_d=p_0+Delta_n;

Blow=4.73e-15 ; %[cm-3s-1] low injection radiative recombination probability
% Values from Altermatt in NUSOD'2005
Bmin=0.2+(0-0.2)/(1+(T/320).^(2.5));
b1=1.5e18+(1e7-1.5e18)/(1+(T/550).^(3.0));
b3=4e18+(1e9-4e18)/(1+(T/365).^(3.54));
Brel=Bmin+(1.00-Bmin)./(1+((n_d+p_d)/b1).^(0.54)+((n_d+p_d)/b3).^(1.25));

B=Brel.*Blow; %[cm-3s-1] radiative recombination probability
Egnarrow=1e-3; % from J. Appl. Phys., Vol. 84, No. 7, 1 October 1998 Andreas Schenk
nieff=ni*exp(Egnarrow/(2*k*T/q));
geeh=1+13*(1-tanh((n_0/3.3e17).^(0.66)));
gehh=1+7.5*(1-tanh((p_0/7e17).^(0.63)));

t_Aug=(Delta_n)./((n_d.*p_d-nieff.^2).*(2.5e-31*geeh*n_0+8.5e-32*gehh*p_0+3e-29*Delta_n.^(0.92)));
t_Rad=(Delta_n)./((n_d.*p_d-nieff.^2).*(B));

end