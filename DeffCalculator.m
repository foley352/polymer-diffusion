function []=DeffCalculator(RH_steps,samplingrate,max_error,data,RH50,RH60,RH70)
%DeffCalculator([RH_steps],[samplingrate],[max_error],[data],[RH50],[RH60],[RH70])
%This matlab function estimates the effective diffusivity of individual RH
%steps by fitting the solution of the differential equation dc/dt = x^-m
%d/dx(x^m Deff*dc/dx) to the experimental data with a boundary condition
%c(edge)=c(Relative Humidity) and a uniform concentration intial
%condition. RH_steps is a vector of humidity steps of the experiment,
%samplingrate is how many time data points to take (1 = all data points, 10 = keep only
%10th data point, etc.).
%RH50, RH60, and RH70 are the humidity step vectors
%for the 50C, 60C and 70C experiments, if different from RH_steps. All inputs
%are optional. Deff is fit be feeding in RH(t) boundary condition data, and
%using the bisection method to get an accurate estimate.
%
%Created by Brandon Foley on 3/21/2024. This work was performed under the
%auspices of the U.S. Department of Energy by Lawrence Livermore National
%Laboratory under Contract DE-AC52-07NA27344.

%check variable inputs, assign default values as necessary
if ~exist('samplingrate')
    samplingrate=1; %how many data points to skip +1; typically 1 (use 10 or 100 if the data file is excessively large to speed up calculations)
end
if ~exist('RH_steps')
    RH_steps=[0:5:90];
end
if ~exist('RH50')
    RH50=RH_steps;
end
if ~exist('RH60')
    RH60=RH_steps;
end
if ~exist('RH70')
    RH70=RH_steps;
end
if ~exist('max_error')
                        max_error=1/2^6;
                    end

%request user for more information
SampleName=input('Sample Name:','s');

%determine geometry of thin lenght scale (i.e., short cylinders are slabs)
geometry=input('Geometry (0=slab, 1=cylinder, 2=sphere):');

%input thickness in cm
thickness=input('Sample Thickness/Radius, cm:');

%if user data is not supplied, the user will be requested to select dynamic
%text files in the form of an IGAsorp or Mettler DVS output file
if ~(exist('data')==1)
    [data,myfiles,path]=parse_files();
end

%convert data to cell if it isn't already
if ~iscell(data)
    data={data};
end

%% extract RH steps from each dynamic file
for j=1:length(data)
    %initialize variables, extract data, find temeprature
    close all
    clear diffusivity ns
    z=data{j};
    z=z(1:samplingrate:end,:);
    T=median(round(data{j}(:,4)));
    
    %determine which RH steps are present in the dynamic file. Could
    %extract these with code but it works best when user-supplied. At
    %higher temperatures (e.g., 50,60, and 70C) the IGAsorp often has
    %slightly different RH (e.g., 6% instead of 5%). This allows the user
    %to specify this difference without having to analyze higher
    %temperatures in a separate batch.
    RH=RH_steps;
    switch T
        case 50
            RH=RH50;
        case 60
            RH=RH60;
        case 70
            RH=RH70;
    end
    clear approxsteps
    
    %Remove nonunique time data points
    [~,x_index]=unique(z(:,1));
    z=z(x_index,:);
    zall=z;
    
    %find where sorption ends
    n2=find(z(:,3)<RH(end)*1.005 & z(:,3)>RH(end)*0.995,1,'last');
    
    %if n2 is empty, try again but don't consider the RH vector input by
    %the user
    if ~min(size(n2))
        n2=find(z(:,3)<max(z(:,3))*1.005 & z(:,3)>max(z(:,3))*0.995,1,'last'); %find where sorption ends
    end
    
    %save the desorption data
    zdesorb=z(max(n2-20,1):end,:);
    
    %crop z so only sorption data is included
    z=z(1:n2,:);
       
    %make RH vector vertical
    RH=reshape(RH,[],1);
   
    %remove outliers using a moving median
    [~,tf]=rmoutliers(z(:,2),'movmedian',3,'ThresholdFactor',2,'SamplePoints',z(:,1));
    z=z(~tf,:);
    
    %define window in which we might reasonably expect noise in RH setpoint
    RHstep=0.005*max(z(:,3)); 
    
    %find where sorption begins
    n1=find(z(:,3)<min(z(:,3)+RHstep),1,'last');
    
    %crop data to include just a few points before the first RH step is reached
    z=z(max(1,n1-2):end,:); 
    
    %subtract out the drymass of the sample
    drymass=min(z(:,2));
    try
        z(:,2)=z(:,2)-drymass;
        zdesorb(:,2)=zdesorb(:,2)-drymass;
        zall(:,2)=zall(:,2)-drymass;
    catch
        warning('Could not initialize drymass')
    end
    
    
    %% Estimate Deff with bisection method
    %initialize variables
    n=1;
    hold on
    counter=1;
    zsorb=z;
    counter1=1;
    counter2=1;
    
    %dir is 0 for sorption, 1 for desorption
    for dir=[0,1]
        if dir
            z=zdesorb;
            RHstore=RH;
            RH=flipud(RH);
        end
        colornum=0;
        for v=1:(length(RH)-1)
            m=RH(v);
            %find when an RH step starts, method depends on
            %sorption/desorption
            if ~dir %sorption
                %find when an RH step begins
                n1=find(z(:,3)<(m+RHstep),1,'last')-5;
                n1=max(n1,1);
            else %desorption
                %find when an RH step begins
                n1=find(z(:,3)>(m-RHstep),1,'last')-5;
                n1=max(n1,1);
            end
            %if we could not find the start of the step, or our subtraction
            %of 5 set the start to negative times, then set n1 = 1. 
            if ~min(size(n1))
                n1=1;
            end
            if n1<1
                n1=1;
            end
            
            %print out RH step
            m=RH(v+1)
            
            %find when an RH step ends, method depends on
            %sorption/desorption
            if ~dir
                n2=find(z(:,3)<m+RHstep & z(:,3)>m-RHstep,1,'last');
            else
                n2=find(z(:,3)>m-RHstep & z(:,3)<m+RHstep,1,'last');
            end
            try
                %if an RH step is wrong sometimes it struggles and tries to skip the step. This corrects this.
                if (n1-ns(n-1))>50
                    n1=ns(n-1);
                    n2=n2-10;
                end
            end
            
            if min(size(n2))
                %store the end index in vector ns
                ns(n)=n2;
                eqbm(n,:)=[m,z(n2,2)];
                %save the final mass at the end of the RH step as the equilibrium point
                n=n+1;
                t=z(n1:n2,1);y=z(n1:n2,2);
                %extract time (t) and mass (y) data for the step
                try
                t=t-t(1);
                catch
                    warning('Time vector is empty. Unable to extract data for RH step')
                end
                RHsSTEP=z(n1:n2,3);
                approxsteps(counter,:)=[z(n1),median(z(n1:n2,3))];
                if approxsteps(counter,2)<0
                    approxsteps(counter,2)=0;
                end
                try
                counter=counter+1;
                t05=(t-t(1)).^0.5; % root time
                t052=sqrt(t05.^2/(abs(median(y(end-5:end))-median(y(1:1)))/drymass));
                
                %normalized mass
                y=(y-median(y(1:1)))/(median(y(end-5:end))-median(y(1:1))); 
                
                %find the points in the data where normalized are in the range of 0.2
                %to 0.5
                n1=find(y<0.2,1,'last')-1;
                n2=find(y<0.5,1,'last')+1;
                if n1<1
                    n1=1;
                end
                if n2>length(y)
                    n2=length(y);
                end
                
                %take a linear fit of these data, and plot. Normalized mass v. t^0.5
                %should be linear in this region if diffusion limited. Slope is used to
                %get a good initial guess for the diffusivity.
                try
                    p=polyfit(t05(n1:n2),y(n1:n2),1);
                    Dh2=p(1)^2*pi()/16;
                    datas=[t,y,RHsSTEP];
                    figure;hold on;plot(t,y,'k','LineWidth',2)
                    yyaxis right;plot(t,RHsSTEP)
                    ylabel('Relative Humidity')
                    yyaxis left
                    xlabel('Time / minutes')
                    ylabel('Relative Mass Change')
                    set(gca,'FontSize',12)
                    box on
                    set(gca,'LineWidth',2)
                    
                    %uses pdepe to fit the solution of the diffusion
                    %equation to the experimental data. Uses bisection
                    %method to find Deff.
                    if Dh2>10/t(end)
                        Dh2=0.1/t(end);
                    end
                    
                    %solve pde for data with RH(t) with real time-dependent
                    %boundary condition for a given D/h^2, output the L2
                    %error between the model and the data
                    %dir0 tells us if our model has too large or too small
                    %of a Dh2
                    [err0,dir0,c]=uptakemodel(Dh2,datas,geometry);
                    plotdata=plot(t,c,'-','LineWidth',2,'Color',lines(1));
                    
                    %Dh20 is the current best guess
                    Dh20=Dh2;
                    
                    %Depending on dir0, we want to shrink or grow our
                    %initial guess
                    if dir0>0
                        pwr=-1;
                    else
                        pwr=1;
                    end
                    
                    %dir2 is the initial positvie or negative residual. We
                    %store that to know when our residual changes sign
                    flags=1;
                    dir2=dir0;
                    iter=0;
                    
                    %while our residual sign remains unchanged, keep
                    %scaling our Dh2 guess by a factor of ten. Once the
                    %residual sign changes we now have an upper and a lower
                    %bound for our best fit. Do not try this more than 3
                    %times.
                    while dir2/dir0>0 & iter<3
                        iter=iter+1;
                        Dh2=Dh2*10^pwr;
                        [err,dir2,c]=uptakemodel(Dh2,datas,geometry);
                        plotdata.XData=t;
                        plotData.YData=c;
                        drawnow()
                    end
                    
                    %if the new bound is better, store that as the new best
                    %guess
                    if err0>err
                       Dh20=Dh2;
                       err0=err;
                    end
                    
                    %keep iterating until the power is less than 1/2^6,
                    %i.e., keep iterating until we have the best fit for D
                    %within a factor of 1/2^6 = 1.56%
                    while abs(pwr)>max_error
                        if dir2>0
                            pwr=-abs(pwr)/2;
                        else
                            pwr=abs(pwr)/2;
                        end
                        Dh2=Dh2*10^pwr;
                        [err,dir2,c]=uptakemodel(Dh2,datas,geometry);
                        
                        %draw the fitting progress
                         plotdata.XData=t;
                         plotdata.YData=c;
                         drawnow()
                         
                        if err0>err
                            Dh20=Dh2;
                            err0=err;
                            plotdata.XData=t;
                            plotdata.YData=c;
                        end
                    end
                    
                    %save picture
                    frame = getframe(gcf);
                    [im, map] = rgb2ind(frame2im(frame),256);
                    imwrite(im, map, ['T_',num2str(T),'_RH',num2str(m),'_',num2str(dir),'.png'], 'png');
                    colornum=colornum+0.05;
                    dataout{counter1}=[t05,y];
                    counter1=counter1+1;               
                    %the RH and slope are stored in a vector
                        diffusivity(n,:)=[m,Dh2];
                catch
                    n=n-1;
                end
                catch
                    warning('Unknown error, probably from saving an image')
                end
            end
            
        end
    end
    RH=RHstore;
    z=[zsorb;zdesorb];
    approxsteps(counter,:)=[z(end,1),approxsteps(end,2)];
    diffusivity=diffusivity(2:end,:);%remove first dummy entry in diffusivity
    
    diffusivity(:,2)=diffusivity(:,2)*(thickness)^2; %calculate effective diffusivity / cm^2 min^-1
   
    diffout{counter2}=[diffusivity(:,1),diffusivity(:,2)];
    
    geos={'Slab','Cylinder','Sphere'};
    %write Deff solutions to files
    if exist('myfiles')
    file=myfiles{j};
    else
        file=sprintf('%s_%0.3f_cm_%s_T_%0.f.txt',SampleName,thickness,geos{geometry+1},T);
    end
    %if there is not a path, save in current directory
    if ~(exist('path')==1)
        path=pwd;
    end
    fid3=fopen(fullfile(path,[file(1:end-4),'_Diffusivity',file(end-3:end)]),'w');
    fprintf(fid3,[SampleName,newline,'Thickness/Radius: ',num2str(thickness),' cm',newline,'Geometry:',geos{geometry+1},newline,'Temperature (',char(176),'C):',num2str(round(median(data{j}(:,4)))),newline,'RH',sprintf('\t'),'D / cm^2/s',sprintf('\t'),'Temperature(',char(176),'C)',sprintf('\t'),'Weights',newline]);
    fprintf(fid3,'%f\t%e\t%f\t%f\n',[diffusivity(:,1)/100,diffusivity(:,2)/60,(round(median(data{j}(:,4)))*ones(length(diffusivity(:,1)),1)),ones(length(diffusivity(:,1)),1)]');
    fclose(fid3);
    
    counter2=counter2+1;
    T=median(z(:,4));
    
    D(j,:)=[median(diffusivity(:,2)),T];

    [~,tf]=min(ns(2:end));
    tf=tf+1;
end
end

%% function for using PDEPE to calculate error, the median residual, and the uptake
function [err,medresid,c]=uptakemodel(Dh2,data,geometry)
    %this functions uses the pdepe solver to estimate the error between the
    %current D/h^2 guess and the experimental data
    %equation is non-dimensionalized so the length scale is x'=x/h or r'=r/R
    %time is non-dimensional as t/(D/h^2) or t/(D/R^2)
xmesh=linspace(0,0.5,101);
if geometry %if cylinder or sphere, xmesh goes from 0 to 1
    xmesh=2*xmesh;
end
%non dimensionalize time
bc(:,1)=data(:,1)*Dh2;

%RH data
bc(:,2)=data(:,3);

%convert column 2 to relative change in RH
bc(:,2)=(bc(:,2)-bc(1,2))/(median(bc(:,2))-bc(1,2));

%make the time step 4x larger than the median time step of the
%experimental data. tstep is equidistant in time which makes interpolating
%the boundary condition significantly faster.
 tstep=median(diff(bc(:,1)))*4;
 
 %new times and boundary conditions for PDEPE
 newxs=[0:tstep:bc(end,1),bc(end,1)]';
 newbc=interp1(bc(:,1),bc(:,2),newxs);
 bc=[newxs,newbc];
 
 %nondimensional time
tspan=data(:,1)*Dh2;

%solve the pde. bc is an array where t' is column 1 and relative change in
%RH is column 2. sol is the c(x,t) data.
sol=pdepe(geometry,@pdefun,@icfun,@(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,bc),xmesh,tspan);

clear bc

%calculate derivative of error with respect to Dh2 (D/h^2), so redo
%everything after multiplying D/h^2 by a 1.01 factor
bc(:,1)=data(:,1)*Dh2*1.01;
bc(:,2)=data(:,3);
bc(:,2)=(bc(:,2)-bc(1,2))/(median(bc(:,2))-bc(1,2));
tstep=median(diff(bc(:,1)))*4;
newxs=[0:tstep:bc(end,1),bc(end,1)]';
newbc=interp1(bc(:,1),bc(:,2),newxs);
bc=[newxs,newbc];
tspan=data(:,1)*Dh2*1.01;

%sol2 is our perturbed c(x,t)
sol2=pdepe(geometry,@pdefun,@icfun,@(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,bc),xmesh,tspan);

%c is the total concentration in our material
c=trapz(xmesh,[sol.*xmesh.^geometry]')/trapz(xmesh,xmesh.^geometry);
%c=sum([(sol.*xmesh.^geometry)]')./sum(xmesh.^geometry);
err=sum((c'-data(:,2)).^2)

c2=trapz(xmesh,[sol2.*xmesh.^geometry]')/trapz(xmesh,xmesh.^geometry);
%c2=sum([(sol2.*xmesh.^geometry)]')./sum(xmesh.^geometry);
err2=sum((c2'-data(:,2)).^2)

%Does increasing Dh2 bring you closer or farther from teh best fit
medresid=(err2>err)*2-1;
end

%% PDE BC and IC functions for quasi-equilibrated Langmuir and Pooling
function [c,f,s]=pdefun(x,t,u,dudx)
    s=0;
    c=1;
    f=dudx;
end

function u0=icfun(x,ic)
u0=0;
end

function [pl,ql,pr,qr]=bcfun(xl,ul,xr,ur,t,bc)
pl=0;
ql=1;
%us actual RH v. time data to set the boundary condition by interpolation
pr=ur-eveninterp1(bc(:,1),bc(:,2),t);
qr=0;
end

%this interpolation function accepts only equally spaced x, which gives
%much fast interpolation than Matlab's interp1 function
function output=eveninterp1(x,y,q)
    %only works for evenly spaced x
    for m=1:length(q)
        n=1+(q(m)-x(1))/(x(2)-x(1));
        a=floor(n);
        b=n-a;
        a=min(a,length(y));
        a=max(a,1);
        if a & ~(a==length(y))
            try
            output(m)=y(a)*(1-b)+y(a+1)*b;
            catch
                1
            end
        elseif ~a
            output(m)=y(1);
        else
            output(m)=y(end);
        end
    end
end

%% parse dynamic files
function [data,myfiles,path]=parse_files()
%select multiple dynamic files in the form of IGAsorp output files
[myfiles,path] = uigetfile('*.txt;*.dat','Select the Input Dynamic File(s)','MultiSelect','on');

%initial file type variable based in dynamic file format
file_type=0;

%store old directory
olddir=pwd;

%check if multiple files selected. Convert string to cell if needed.
if ~iscell(myfiles)
    myfiles={myfiles};
end
cd(path)
mkdir figures

%evaluates which version of IGAsorp file is input, and adjusts the column
%number and number of lines to skip as necessary. Not guaranteed to work
%with all IGAsorp files. Some adjustments may be required to validate all
%data inputs. Ultimately, the final results is an array, "data", that has
%4 columns: [Time (s), Mass (mg), RH (%), Temperature (C)]. This code can
%be replaced with anything that extracts the data to this format.
for  j=1:length(myfiles)
    file=myfiles{j};
    fid1=fopen(file);
    if fid1>-1
        flag=0;
        while flag==0 & ~feof(fid1)
            tline=fgetl(fid1);
            %does the file contain "Time"
            IndexC = strfind(tline,'Time');
            if ~min(size(IndexC))
                %does the file contain "TIME"
                IndexC = strfind(tline,'TIME');
                if min(size(IndexC))
                    %Does the file containt a column number "A3"
                    IndexC = strfind(tline,'A3');
                    if min(size(IndexC))
                        %then the number of columns in this file is 11
                        columns=11;
                        %Does the file contain a column number "A2"
                    elseif min(size(strfind(tline,'A2')))
                        %then the number of columns in this file is 10 and
                        %the file_type is 1
                        columns=10;
                        file_type=1;
                    else
                        %if "A2" and "A3" are not found in these files, the
                        %the column number is 7 and the file_type is 1
                        columns=7;
                        file_type=1;
                    end
                    %skip a line and exit while loop
                        tline=fgetl(fid1);
                        flag=1;
                end
            else
                %check for "A =" line, if present file is an IGAsorp file,
                %if absent, it is a Mettler DVS file
                IndexC = strfind(tline,'A =');
                if min(size(IndexC))
                    %this IGAsorp file has 25 columns. Skip 26 lines to get
                    %to data array in file.
                    flag=1;
                    columns=25;
                    for v=1:26
                        tline=fgetl(fid1);
                    end
                else
                    %check for "in hours" in header, then it is a Mettler
                    %DVS file
                    IndexC = strfind(tline,'in hours');
                    if min(size(IndexC))
                        %Mettler DVS files have 7 columns, skip 1 line to
                        %get to data array
                        flag=1;
                        columns=7;
                        tline=fgetl(fid1);
                    end
                end
            end
        end
        %if you reach the end of the file, throw an error to user
        if feof(fid1)
            error("Could not parse dynamic file type")
        end
        %extract data from file with columns specified above
        data{j}=fscanf(fid1,'%f',[columns,inf])';
        
 %% extract data
 %define data as Time (minutes), Mass (mg), Relative Humidity (%),
 %Temperature (C)
        %if file has 11 columns or is file_type 1 then extract the data in
        %this format.
        if columns==11 | file_type %historical IGAsorp file
            data{j}=data{j}(:,[1:4]);
            if exist('drymass')
                data{j}(:,2)=data{j}(:,2)+drymass;
            end
            %if the file has 25 columns, we want to extract columns 1-3 and
            %14
        elseif columns==25 %Modern IGAsorp file
            data{j}=data{j}(:,[1:3,14]);
            %if the file has 7 columns, we want columns, 1,3,4, and 7
        elseif columns==7 %Mettler DVS
            data{j}=data{j}(:,[1,3,4,7]);
            %convert hours to minutes
            data{j}(:,1)=data{j}(:,1)*60;
        else
            error('Unclear file format')
        end
        fclose(fid1)
    end
end
end