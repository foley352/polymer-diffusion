function [iso_data,Deff_data]=FitIsoAndDeff()
%[iso_data,Deff_data]=FitIsoAndDeff()
%Find best fit of isotherm and Deff data simulataneously for
%sorption-diffusion of species through materials. The isotherms are assumed
%to be the sum of mobile Langmuir, mobile Henry's mobile pooling, immobile
%Langmuir, immobile Henry's, immobile pooling:
%c = bm*a/(1+bm*a)*Lm + Hm*a + pm*a^nm + bi*a/(1+bi*a)*Li + Hi*a + pi*a^ni 
%where parameters "m" are mobile and "i" are immobile. a = activity
%it is assume that only mobile species diffuse and that
%non-linearities in mobile concentration v. activity imply non-constant
%activity coefficient. Deff is thus approximated by the function
%Deff ~ c_m/a/(d(c_m+c_i)/da)*D
%where c_m(a) is a function of mobile concentration versus activity, c_i(a)
%is the immobile concentration as a function of a, and D is a function only
%of temperature
%this matlab function requests users to (1) specify mobile and immobile
%desorption modes and (2) select data files containing all isotherms and
%all effective diffusivities so that multiple temperatures can be fit
%simultanouesly at the same time. fminunc is used to find optimal fits for
%100 initial guesses generated by the latin-hypercube method. The best fit
%of these is considered the global maximum. Confidence t-intervals are
%calcualted and the best fit is output to a .txt file. van't Hoff and
%Arrhenius temperature dependencies are assumed for equilibrium and
%diffusivity parameters, respectively.
%Effective diffusivity input files must at a minimum contain a header
%containt D / cm^2/s and one row of data, for example:
%
%RH	D / cm^2/s	Temperature(°C)	Weights
%0.050000	3.017512e-05	30.000000	1.000000
%
%where Weights is the weighting factor for this data point in the fit (set
%outliers to 0 for improved fits)
%
%Isotherm input files must include "Sorption" and "Desorption" Headers and
%formatted as below for example (additional header lines are ok):
%
% Sorption
% RH (%)	Mass (mg)	Uptake (mg/g)	Temperature (°C)
% 0.000000	953.798200	0.000000	50.000000
% 0.000000	953.799100	0.000944	50.000000
% 5.018800	954.224700	0.447160	50.000000
% 9.989300	954.393300	0.623927	50.000000
% 15.011500	954.525000	0.762006	50.000000
% 19.980000	954.641200	0.883835	50.000000
% 24.997800	954.746400	0.994131	50.000000
% 30.021000	954.848800	1.101491	50.000000
% 34.994500	954.941500	1.198681	50.000000
% 40.014200	955.027200	1.288533	50.000000
% 44.987700	955.111100	1.376497	50.000000
% 50.007400	955.186300	1.455339	50.000000
% 54.977500	955.260700	1.533343	50.000000
% 59.993700	955.327800	1.603694	50.000000
% 65.023300	955.396200	1.675407	50.000000
% 69.988400	955.463400	1.745862	50.000000
% 75.012100	955.528700	1.814325	50.000000
% 79.984100	955.597700	1.886667	50.000000
% 85.006800	955.660700	1.952719	50.000000
% 89.976300	955.725700	2.020868	50.000000
% 
% Desorption
% RH (%)	Mass (mg)	Uptake (mg/g)	Temperature (°C)
% 89.976300	955.725700	2.020868	50.000000
% 89.978800	955.734300	2.029884	50.000000
% 84.999400	955.669200	1.961631	50.000000
% 79.985100	955.601900	1.891071	50.000000
% 75.013100	955.540500	1.826697	50.000000
% 69.988900	955.476600	1.759701	50.000000
% 65.022400	955.412100	1.692077	50.000000
% 59.997700	955.345500	1.622251	50.000000
% 54.975500	955.278900	1.552425	50.000000
% 50.008900	955.202900	1.472743	50.000000
% 44.986200	955.127000	1.393167	50.000000
% 40.014200	955.045900	1.308138	50.000000
% 34.990500	954.957700	1.215666	50.000000
% 30.021500	954.863600	1.117008	50.000000
% 25.001700	954.766500	1.015204	50.000000
% 19.980000	954.658300	0.901763	50.000000
% 15.012500	954.546200	0.784233	50.000000
% 9.988300	954.408300	0.639653	50.000000
% 5.017300	954.240700	0.463935	50.000000
% 0.000000	953.877700	0.083351	50.000000
%
%Created by Brandon Foley on 3/21/2024. This work was performed under the
%auspices of the U.S. Department of Energy by Lawrence Livermore National
%Laboratory under Contract DE-AC52-07NA27344.


%request user inputs to get sample name, mobiel and immobile sorption modes
SampleName=input('Sample Name:','s');
modesm=input('MOBILE sorption modes: Langmuir (1/0) Henrys (1/0) Pooling (1/0) entered as [L,H,P]:');
modesi=input('IMMOBILE sorption modes: Langmuir (1/0) Henrys (1/0) Pooling (1/0) entered as [L,H,P]:');
TrefK=293.15;

%initialize variables
TrefC=TrefK-273.15;
bmin=0.1;
nmin=1;
lf=modesm(1); %0 or 1 to include langmuir
hf=modesm(2); %0 or 1 to include henry's
pf=modesm(3); %0 or 1 to include pooling
lfi=modesi(1); %0 or 1 to include langmuir
hfi=modesi(2); %0 or 1 to include henry's
pfi=modesi(3); %0 or 1 to include pooling

%not used, but could be modified to allow Lanmguir sorption mode to have a
%linear temperature tempendence
dldtf=0; %0 or 1 to include temperature dependent langmuir

%user input for selecting multiple iso and diffusivitiy files
[myfiles,path] = uigetfile('*.txt;*.dat','Select the Input Iso File(s)','MultiSelect','on'); %select multiple files that you want plots of; all must either be SD or iso files
[myfiles2,path2] = uigetfile('*.txt;*.dat','Select the Input Diffusivity File(s)','MultiSelect','on'); %select multiple files that you want plots of; all must either be SD or iso files

%convert myfiles to cell if only one file was selected
if ~iscell(myfiles)
    myfiles={myfiles};
end
if ~iscell(myfiles2) & myfiles2
    myfiles2={myfiles2};
end

%save old directory and change directory to path of iso files
olddir=pwd;
cd(path)

%iterate through seach iso file, extract data, and store in cell "iso_data"
for j=1:length(myfiles)
    file=myfiles{j};
    fid1=fopen(file)
    flag=0;
    while flag==0 & ~feof(fid1)
        tline=fgetl(fid1);
        IndexC = strfind(tline,'Sorption');
        if min(size(IndexC))
            tline=fgetl(fid1);
            iso_data{j}=fscanf(fid1,'%f',[4,inf])';
        end
        IndexC = strfind(tline,'Desorption');
        if min(size(IndexC))
            tline=fgetl(fid1);
            A=fscanf(fid1,'%f',[4,inf])';
            iso_data{j}=[iso_data{j};A];
            flag=1;
        end
    end
    fclose(fid1)
end

%iterate through dynamic files, extract diffusivity, store in cell "Deff_data"
if iscell(myfiles2)
    for j=1:length(myfiles2)
        file=myfiles2{j};
        fid1=fopen(fullfile(path2,file),'r');
        flag=1;
        while ~feof(fid1) & flag
            tline=fgetl(fid1);
            if min(size(strfind(tline,'D / cm^2/s')))
                flag=0;
            end
        end
        A=fscanf(fid1,'%f',[4,inf]);
        Deff_data{j}=A';
    end
end

%variable storing which modes are mobile and immobile--details for fitting
%functions
ns=[lf,lf,hf,pf,pf,lf,hf,pf,lfi,lfi,hfi,pfi,pfi,lfi,hfi,pfi,1,1];

%generate 100 latinhypercube intial guesses
guesses=latinhypercube(100,sum(ns));
guesses(:,~~ns)=guesses;

%converts random numbers to roughly constrained bounds for each initial
%guess. These variables are actually the natural log of the parameters
%these guesses include mobile and immobile b, L, H, p, n and Ea,B Ea,H, and
%Ea,p, D, and Ea,D.
guesses(:,1)=(guesses(:,1)-0.5)*10;
guesses(:,2)=(guesses(:,2)-0.5)*20;
guesses(:,3)=(guesses(:,3)-0.5)*10;
guesses(:,4)=(guesses(:,4)-0.5)*10;
guesses(:,5)=guesses(:,5)*exp(1);
guesses(:,[6:8])=(guesses(:,[6:8])-0.5)*10;

guesses(:,9)=(guesses(:,9)-0.5)*10;
guesses(:,10)=(guesses(:,10)-0.5)*20;
guesses(:,11)=(guesses(:,11)-0.5)*10;
guesses(:,12)=(guesses(:,12)-0.5)*10;
guesses(:,13)=(guesses(:,13))*exp(1);
guesses(:,[14:16])=(guesses(:,[14:16])-0.5)*10;

guesses(:,17)=(guesses(:,17)*1.5+0.5);
guesses(:,18)=(guesses(:,18)*1.5+0.5);

%fit Deff data assuming Deff has an arrhenius dependence to get a really
%good intial guess for Deff arrhenius parameters, use these paramters as
%the initial guess for each optimization
for j=1:length(Deff_data)
    Dvec(j)=median(Deff_data{j}(:,2));
    Tvec(j)=median(Deff_data{j}(:,3));
end

p=polyfit(1./(Tvec+273.15),log(Dvec),1);
guesses(:,17)=polyval(p,1/TrefK)*(guesses(:,17));;
guesses(:,18)=p(1)/1000*(guesses(:,18));
guesses(:,~ns)=[];

%initial variables and set optimizaiton options
besterr=1e10;
options = optimoptions('fminunc','Display','off','MaxFunctionEvaluations', inf,'MaxIterations',1000);
RH=0;Ts=0;ys=0;
for j=1:length(iso_data)
    RH=[RH;iso_data{j}(:,1)];
    Ts=[Ts;iso_data{j}(:,4)];
    ys=[ys;iso_data{j}(:,3)];
end
RH=RH/100;
RH(1)=[];
Ts(1)=[];
ys(1)=[];

%cycle through all intiial guesses and minimize objective function, store
%the best fits in xsave and lowest error as besterr
tic
for iter=1:length(guesses)
    try
    xguess=guesses(iter,:);
    nsum=cumsum(ns);
    [x,b,exitflag]=fminunc(@(x) abs(triplemodefitanddiff([exp(x(max(1,nsum(1:5))))+[bmin,0,0,0,nmin],x(max(1,nsum(6:8))),exp(x(max(1,nsum(9:13))))+[bmin,0,0,0,nmin],x(max(1,nsum(14:16))),exp(x(nsum(17))),x(nsum(18))],lf,hf,pf,lfi,hfi,pfi,dldtf,iso_data,Deff_data,TrefK)),xguess,options);
    if exitflag<1
        warning('Did not finish minimizing for THIS particular guess')
    end
    if b<besterr
        besterr=b;
        sprintf('Best Iteration:%0.f, Lowest Error:%f',[iter,besterr])
        xsave=x;
    end
    catch
        warning("Encountered error when trying to minimize this initial guess")
    end
end

%store best fit as x. x(~~ns)=x is here becuase fminunc has a variable
%number of parameters depending on how many mobile/immobile modes the user
%included
x=xsave;
x(~~ns)=x;
%convert x (logarithmic variables with lower bounds) to actual variabeles
%x2
x2=[exp(x(1:5))+[bmin,0,0,0,nmin],x(6:8),exp(x(9:13))+[bmin,0,0,0,nmin],x(14:16),exp(x(17)),x(18)].*[lf,lf,hf,pf,pf,lf,hf,pf,lfi,lfi,hfi,pfi,pfi,lfi,hfi,pfi,1,1];

%if one of the varaiable fits is very large at the upperbound of where
%matlab converts floating point numbers to infinity (realmax), then
%multiply this very large parameter by 0.9 to avoid "inf" errors. This
%likely isn't a great fit if this is required
x2=x2.*(0.9.^isinf(x2*1.001));
    
%calculate baseline error and model fits to iso (fit11) and Deff (fit12)
[err,fit11,fit12]=triplemodefitanddiff(x2,lf,hf,pf,lfi,hfi,pfi,dldtf,iso_data,Deff_data,TrefK);
dataout2=fit11;dataout3=fit12;

%use function parsefit2 to combine our model fits into a single vector
fit1=parsefit2(fit11,fit12);

%iterate through each parameter and purturb by 0.1% to calculate the
%sensitivity of the fit to that parameter. This will be used to caluclate
%the Hessian matrix and estimate the confidence t-intervals
for j=1:length(x2)
    dx=0.001*x2(j);
    x2(j)=1.001*x2(j);
    [~,fit21,fit22]=triplemodefitanddiff(x2,lf,hf,pf,lfi,hfi,pfi,dldtf,iso_data,Deff_data,TrefK);
    fit2=parsefit2(fit21,fit22);
    X2(:,j)=(fit2-fit1)/(dx);
    x2(j)=x2(j)/1.001;
end

%for historical reasons, the triplemodefitanddiff function multiplied
%activation energies by 1000, so the actual parameters are 1000x larger
%than the fit gives and the sensitivities (X2) are 1000x less
x2([6:8,14:16,18])=x2([6:8,14:16,18])*1000;
X2(:,[6:8,14:16,18])=X2(:,[6:8,14:16,18])/1000;

%x2 and X2 contain all possible mobile and immobile parameters, but some of
%these are set to zero by the user. This code removes theses
ns=~[lf,lf,hf,pf,pf,lf,hf,pf,lfi,lfi,hfi,pfi,pfi,lfi,hfi,pfi,1,1];
xkeep2=x2.*[~ns];
x2=x2(~ns);

%this is the Jacobian
X2=X2(:,~ns);

%estimate the inverse of the Hessian
M=inv(X2'*X2);

%number of parameters
params=length(x2);

%degrees of freedom
dof=length(X2)-params;

%standard error squared
se2=(err/dof);

%t-distributions
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));                                % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;                                          % 1-tailed t-distribution
t_inv = @(alpha,v) fzero(@(tval) (max(alpha,(1-alpha)) - tdist1T(tval,v)), 5);  % T-Statistic Given Probability ‘alpha’ & Degrees-Of-Freedom ‘v’

%tvalue of 95% confindence interval for 2-tailed distribution with dof
%degrees of freedom
tval=t_inv(0.05/2,dof);

%xerr2 is the error in each parameter
xerr2=tval*sqrt(se2*diag(M));

%xci2 is the confidence interval for each parameter
xci2=[x2'-xerr2,x2'+xerr2];

%if user included Deff datasets, make some figures showing data fits
if exist('Deff_data')
m=1;
xerr22=xerr2;
for j=1:length(xkeep2)
    if xkeep2(j)
        xcis2(j,:)=xci2(m,:);
        xerr2(j,:)=xerr22(m,:);
        m=m+1;
    else
        xcis2(j,:)=[0,0];
        xerr2(j,:)=[0];
    end
end
figure;hold on

legendstr={};
for j=1:length(dataout3)
    colors=lines(length(dataout3));
    plot(([0;dataout3{j}(1:end-1,1)]+dataout3{j}(1:end,1))/2,dataout3{j}(:,2),'o','MarkerFaceColor',colors(j,:),'Color',colors(j,:));
    plot(([0;dataout3{j}(1:end-1,1)]+dataout3{j}(1:end,1))/2,dataout3{j}(:,5),'-','Color',colors(j,:),'LineWidth',2);
    legendstr{2*j-1}=['T = ',num2str(median(dataout3{j}(:,3))),' ',char(176),'C'];
    legendstr{2*j}=['T = ',num2str(median(dataout3{j}(:,3))),' ',char(176),'C model'];
end
box on
set(gca,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('Relative Humidity')
ylabel('D_{eff} / cm^2 s^{-1}')
legend(legendstr,'Location','Best')
saveas(gcf,[SampleName,'_DiffusivityFits'],'fig');
saveas(gcf,[SampleName,'_DiffusivityFits'],'png');
end

%if user included Deff datasets, store results in files
if exist('Deff_data')
fid1=fopen('Iso+Dyn Parameters.txt','w')
fprintf(fid1,['uptake / mg g^-1 = b*RH/(1+b*RH)*L + H*RH + a*RH^n for mobile (m) and immobile (i) species',newline]);
fprintf(fid1,['bm(T=',num2str(TrefC),' ',char(176),'C)',sprintf('\t'),'Lm',sprintf('\t'),'Hm(T=',num2str(TrefC),' ',char(176),'C)',sprintf('\t'),'am(T=',num2str(TrefC),' ',char(176),'C)',sprintf('\t'),'n',sprintf('\t'),'dlnbm/d1/T',sprintf('\t'),'dlnHm/d1/T',sprintf('\t'),'dlnam/d1/T',sprintf('\t'),'bi(T=',num2str(TrefC),' ',char(176),'C)',sprintf('\t'),'Li',sprintf('\t'),'Hi(T=',num2str(TrefC),' ',char(176),'C)',sprintf('\t'),'ai(T=',num2str(TrefC),' ',char(176),'C)',sprintf('\t'),'n',sprintf('\t'),'dlnbi/d1/T',sprintf('\t'),'dlnHi/d1/T',sprintf('\t'),'dlnai/d1/T',sprintf('\t'),'D(T=',num2str(TrefC),' ',char(176),'C)',sprintf('\t'),'dlnD/d1/T',newline]);
fprintf(fid1,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%f\t',xkeep2');
fprintf(fid1,'\n\n');
fprintf(fid1,['95%% confidence interval\n']);

fprintf(fid1,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%f\t\n',xcis2);
fprintf(fid1,'\n\n');
fprintf(fid1,['95%% error\n']);
fprintf(fid1,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%f\t\n',xerr2);

fprintf(fid1,'\n\n');
fprintf(fid1,['95%% relative error\n']);
fprintf(fid1,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n',abs(xerr2./xkeep2'));

fclose(fid1);
fclose all
x=xkeep2;
else
    error('Must include diffusivity data')
end

%make parity plots of isotherm fits
figure;hold on;colors=lines(length(iso_data));
parity=[0,0];

b0m=x(1);
L0m=x(2);
H0m=x(3);
a0m=x(4);
n0m=x(5);
dlnbdinvTm=x(6);
dlnHdinvTm=x(7);
dlnadinvTm=x(8);

b0i=x(9);
L0i=x(10);
H0i=x(11);
a0i=x(12);
n0i=x(13);
dlnbdinvTi=x(14);
dlnHdinvTi=x(15);
dlnadinvTi=x(16);

D0=x(17);
dlnDdinvT=x(18);

dLdT=0;
dndT=0;

for j=1:length(iso_data)
    RHs=iso_data{j}(:,1)/100;
    T=median(iso_data{j}(:,4));
    legendstrings{2*j-1}=['T = ',num2str(round(T/10)*10),' ',char(176),'C'];
    legendstrings{2*j}=['Fit'];
    
    bm=b0m*exp(dlnbdinvTm*(1/(T+273.15)-1/TrefK));
    Lm=L0m+dLdT*(T-TrefC);
    Hm=H0m*exp(dlnHdinvTm*(1/(T+273.15)-1/TrefK));
    am=a0m*exp(dlnadinvTm*(1/(T+273.15)-1/TrefK));
    nm=n0m+dndT*(T-TrefC);
    
    bi=b0i*exp(dlnbdinvTi*(1/(T+273.15)-1/TrefK));
    Li=L0i+dLdT*(T-TrefC);
    Hi=H0i*exp(dlnHdinvTi*(1/(T+273.15)-1/TrefK));
    ai=a0i*exp(dlnadinvTi*(1/(T+273.15)-1/TrefK));
    ni=n0i+dndT*(T-TrefC);
    
    iso_data{j}(:,5)=bm*RHs./(1+bm*RHs)*Lm*lf+Hm*RHs*hf+pf*am*RHs.^nm+bi*RHs./(1+bi*RHs)*Li*lfi+Hi*RHs*hfi+pfi*ai*RHs.^ni;
    plot(iso_data{j}(:,1)/100,iso_data{j}(:,3),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:))
    fplot(@(RHs) bm*RHs./(1+bm*RHs)*Lm*lf+Hm*RHs*hf+pf*am*RHs.^nm+bi*RHs./(1+bi*RHs)*Li*lfi+Hi*RHs*hfi+pfi*ai*RHs.^ni,[0,max(RHs)],'Color',colors(j,:),'LineWidth',2);
    parity=[parity;iso_data{j}(:,3),iso_data{j}(:,5)];
end

parity=parity(2:end,:);

%calculate rsquare and adjusted rsquare
rsquare=1-sum((parity(:,2)-parity(:,1)).^2)./sum((parity(:,1)-mean(parity(:,1))).^2);
rsqradj=1-((1-rsquare)*(length(parity)-1)/(length(parity)-lf*3-hf*2-pf*3-1));

%format figure and save
nx=ylim;
box on
xlabel('Relative Humidity')
ylabel('Uptake / mg g^{-1}')
legend(legendstrings,'Location','Northwest','NumColumns',2)
ny=ylim;ny(1)=0;
nx=xlim;nx(1)=0;
xlim(nx);
ylim(ny);
saveas(gcf,'IsothermFits','fig');
saveas(gcf,'IsothermFits','png');
figure;hold on
for j=1:length(iso_data)
    plot(iso_data{j}(:,5),iso_data{j}(:,3),'o','MarkerFaceColor',colors(j,:))
end
nx=xlim;ny=ylim;
hold on;fplot(@(x) x,nx,'LineWidth',2,'MarkerFaceColor',colors(1,:));
xlim([0,nx(2)]);
ylim([0,ny(2)]);
box on

set(gca,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('Predicted Uptake / mg g^{-1}')
ylabel('Measured Uptake / mg g^{-1}')
saveas(gcf,'ParityPlot','fig');
saveas(gcf,'ParityPlot','png');

cd(olddir)
end

function [err,data,data2]=triplemodefitanddiff(x,lf,hf,pf,lfi,hfi,pfi,dldtf,data,data2,TrefK,thickness)
TrefC=TrefK-273.15;
err=0;
%mobile parmeters
b0m=x(1);
L0m=x(2);
H0m=x(3);
a0m=x(4);
n0m=x(5);

%dlnbinvT is d ln b / d(1/T) = Ea,B/R, multiplied here by 1000 to scale the
%values to be nice numbers. This is corrected in fit results later

dlnbdinvTm=x(6)*1000;
dlnHdinvTm=x(7)*1000;
dlnadinvTm=x(8)*1000;

%immobile parameters
b0i=x(9);
L0i=x(10);
H0i=x(11);
a0i=x(12);
n0i=x(13);

%dlnbinvT is d ln b / d(1/T) = Ea,B/R, multiplied here by 1000 to scale the
%values to be nice numbers. This is corrected in fit results later
dlnbdinvTi=x(14)*1000;
dlnHdinvTi=x(15)*1000;
dlnadinvTi=x(16)*1000;

D0=x(17);
dlnDdinvT=x(18)*1000;

dLdT=0;
dndT=0;

%iterate through each temperture of isotherm dataset and calculate square
%error of model to experiment
for j=1:length(data)
    try
    RHs=data{j}(:,1)/100;
    T=median(data{j}(:,4));
    
    bm=b0m*exp(dlnbdinvTm*(1/(T+273.15)-1/TrefK));
    Lm=L0m+dLdT*(T-TrefC);
    Hm=H0m*exp(dlnHdinvTm*(1/(T+273.15)-1/TrefK));
    am=a0m*exp(dlnadinvTm*(1/(T+273.15)-1/TrefK));
    nm=n0m+dndT*(T-TrefC);
    
    bi=b0i*exp(dlnbdinvTi*(1/(T+273.15)-1/TrefK));
    Li=L0i+dLdT*(T-TrefC);
    Hi=H0i*exp(dlnHdinvTi*(1/(T+273.15)-1/TrefK));
    ai=a0i*exp(dlnadinvTi*(1/(T+273.15)-1/TrefK));
    ni=n0i+dndT*(T-TrefC);
    
    data{j}(:,5)=bm*RHs./(1+bm*RHs)*Lm*lf+Hm*RHs*hf+pf*am*RHs.^nm+bi*RHs./(1+bi*RHs)*Li*lfi+Hi*RHs*hfi+pfi*ai*RHs.^ni;
    err=err+1e2*sum((data{j}(:,5)-data{j}(:,3)).^2)/median(data{j}(:,3))^2;
    catch
        warning("error calculating model isotherm outputs with these parameters (continuing on to new initial guess)")
    end
end

%iterate through each Deff dataset, calculate errors
for j=1:length(data2)
    try
    RHs=[0;data2{j}(:,1)];
    T=median(data2{j}(:,3));
    
    bm=b0m*exp(dlnbdinvTm*(1/(T+273.15)-1/TrefK));
    Lm=L0m+dLdT*(T-TrefC);
    Hm=H0m*exp(dlnHdinvTm*(1/(T+273.15)-1/TrefK));
    am=a0m*exp(dlnadinvTm*(1/(T+273.15)-1/TrefK));
    nm=n0m+dndT*(T-TrefC);
    
    bi=b0i*exp(dlnbdinvTi*(1/(T+273.15)-1/TrefK));
    Li=L0i+dLdT*(T-TrefC);
    Hi=H0i*exp(dlnHdinvTi*(1/(T+273.15)-1/TrefK));
    ai=a0i*exp(dlnadinvTi*(1/(T+273.15)-1/TrefK));
    ni=n0i+dndT*(T-TrefC);

    D=D0*exp(dlnDdinvT*(1/(T+273.15)-1/TrefK));
    
    RHs=0.5*(RHs(1:end-1)+RHs(2:end));
    mobile=bm*RHs./(1+bm*RHs)*Lm*lf+Hm*RHs*hf+pf*am*RHs.^nm;
    dmdRH=bm*RHs./(1+bm*RHs).^2*Lm*lf+Hm*RHs*hf+pf*am*nm*RHs.^nm;
    didRH=bi*RHs./(1+bi*RHs).^2*Li*lfi+Hi*RHs*hfi+pfi*ai*ni*RHs.^ni;
    
    

    data2{j}(:,5)=D.*mobile./(dmdRH+didRH);
    try
        err=err+sum((data2{j}(:,2)-data2{j}(:,5)).^2.*data2{j}(:,4))/median(data2{j}(:,2))^2*1e0;
    catch
        warning("error calculating model diffusivity outputs with these parameters (continuing on to new initial guess)")
    end
    catch
        warning("error calculating model diffusivity outputs with these parameters (continuing on to new initial guess)")
    end
end
end

%generates n intiial guess for m parameters using the latin hypercube
%method. varibles scale between 0 and 1
function x=latinhypercube(n,m)
%n = number of row, m = number of columns in a latinhypercube generated
%matrix
x=zeros(n,m);
for q=1:n
    x(q,:)=rand(1,m)/n+(q-1)/n;
end
for q=1:m
    [~,idx]=sort(rand(n,1));
    x(:,q)=x(idx,q);
end
end

%reorganizes isotherm and Deff fits into a usable format for later
%calcuations
function output=parsefit2(fit1,fit2)
output=0;
for j=1:length(fit1)
    output=[output;fit1{j}(:,5)/median(fit1{j}(:,3))*sqrt(1e2)];
end
for j=1:length(fit2)
    output=[output;fit2{j}(:,5)/median(fit2{j}(:,2)).*sqrt(1e0*fit2{j}(:,4))];
end
output=output(2:end);
end