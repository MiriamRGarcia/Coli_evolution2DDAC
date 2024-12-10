clear all
close all

%%%=========================================================================
%%%                 Load data
%%%=========================================================================

%% load data
load data_vial1


%% remove last half day (probably with biofilms)

ind_r   = find(Timeh>=60);                                                  %h
t       = Timeh(1:ind_r(1)-1);                                              %h

v.uFm   = Mediapumponoff(1:ind_r(1)-1);                                     % 0,1
v.uFd   = Drugpumponoff(1:ind_r(1)-1);                                      % 0,1
ind_uFm = find(v.uFm==1);                       
ind_uFd = find(v.uFd==1);

OD      = OD(1:ind_r(1)-1);                                                 % OD
Ca      = DruginvialmgL(1:ind_r(1)-1);                                      % ppm


%% filter sampling points (too often
p.dilt0	= 1.76;                                                             % time starting adding medium
inddilu = find(t>p.dilt0);
dt      = 12;                                                               % 12 1min 2*12*2min
ind=[inddilu(1) 1:60:size(v.uFm,1)  ind_uFd' ind_uFd'+dt ind_uFm']';

ind=unique(ind);
ind=sort(ind);
t   = t(ind);                                                               %h
v.uFm = v.uFm(ind);
v.uFd = v.uFd(ind);
OD  = OD(ind);
Ca  = Ca(ind);

%% parameters
p.V     = 12;                                                               % L
p.Cin   = 40;                                                               % mg/L
p.Fm    = (1/5)*60*60;                                                      % 1mL per minute in paper, but now in 5 sg 1ml
p.Fd    = (1/5)*60*60;                                                      % 1mL per minute in paper,(1mL/1min/60) but to add 1 mL in 5 seconds (used for measurements... => 1mL 1/(5/60)min/60.... 60mL en una hora, en 5sg 0.083333
p.addedml= 1;                                                         % mL added when pump is activated
init    = 1e8;
nr      = 3;

est= [1.7741       2.7457     0.023957      0.20435       0.2264       1.8202            3   9.7258e-05       3.2196       1.5897      0.57444      0.36199   1.0017e-05       1.0639       9.0363];


p.G     = est(1);
p.ks    = est(2);
p.YY    = est(3); %
p.D     = est(4);
p.n     = est(5);
p.G2	= est(6);
p.ks2	= est(7);
p.D2	= est(8);
p.n2	= est(9);
p.G3	= est(10);
p.ks3	= est(11);
p.D3	= est(12);
p.n3	= est(13);
p.mu1	= 10^-est(14);
p.mu2	= 10^-est(15);




%% init condition and pumps initialization

X0=[init 0 0 1 0 ];

YY(1,:)=X0;
ii=find(v.uFm==1);
v.uFm(ii(1))=0;
ii=find(v.uFd==1);
v.uFd(ii(1))=0;

%% loop for integrating ODEs
for ii=1:size(t,1)-1
    if v.uFd(ii)>0
        Vtod=p.V+v.uFd(ii)*p.addedml;
        Vadd=v.uFd(ii)*p.addedml;
        X0(1:nr) =(p.V*X0(1:nr)+0*Vadd)./Vtod;
        X0(end-1)=(p.V*X0(end-1)+1*Vadd)/Vtod;
        X0(end)  =(p.V*X0(end)+40*Vadd)/Vtod;
    end

    if v.uFm(ii)>0
        Vtom=p.V+v.uFm(ii)*p.addedml;
        Vadm=v.uFm(ii)*p.addedml;
        X0(1:nr) =(p.V*X0(1:nr)+0*Vadm)./Vtom;
        X0(end-1)=(p.V*X0(end-1)+1*Vadm)/Vtom;
        X0(end)  =(p.V*X0(end)+0*Vadm)/Vtom;
    end
    
    [ti,yi]=ode15s(@theODEmodel,[t(ii) t(ii)+(t(ii+1)-t(ii))/2 t(ii+1)],X0,[],p);
    
    YY(ii+1,:)=yi(end,:);
    X0=yi(end,:);
end

XX=YY(:,1:nr);
NN=sum(XX,2);
S=YY(:,nr+1);
C=YY(:,nr+2);
bet=(0.012)/init;
ODm=bet*NN;



%% plot results
figure
set(gcf,'OuterPosition',[998         305        1055        1113])

subplot(3,1,1),
plot(t,OD,'g');hold on
plot(t,ODm,'k')
legend('Measured OD','Model Simulation of OD')
ylabel('OD')
title('Model vs Data')

subplot(3,1,2),
plot(t,sum(XX,2),'k');hold on
plot(t,XX(:,1),	'Color',[1 1 0]),hold on
plot(t,XX(:,2),'Color',[1 0.5 0]),hold on
plot(t,XX(:,3),'Color',[1 0 0]),hold on
legend('Total cells','Susceptible cells','Intermediate cells','Resistant cells')
ylabel('cfu/mL')
title('Dynamics of subpopulations')

subplot(3,1,3),
[ax,h1,h2] =plotyy(t,C,t,S);
set(h1, 'Color', 'm');
set(h2, 'Color', 'b');
ylabel(ax(1), 'DDAC [ppm]');
ylabel(ax(2), 'Enviromental stress');
xlabel('Time (h)')
title('Stresse')


function dx=theODEmodel(t,x,p)

X=x(1);
X2=x(2);
X3=x(3);
S=x(end-1);
C=x(end);

dX=p.G*X*S^p.ks-p.D*X*C^p.n-p.mu1*X;%;
dX2=p.G2*X2*S^p.ks2-p.D2*X2*C^p.n2+p.mu1*X-p.mu2*X2;%;
dX3=p.G3*X3*S^p.ks3-p.D3*X3*C^p.n3+p.mu2*X2;%;
dS=-(p.YY/1e8)*(p.G*X*S^p.ks+p.G2*X2*S^p.ks2+p.G3*X3*S^p.ks3);
dC=0;

dx=[dX;dX2;dX3;dS;dC];    

end





