close all ;
clear all ;



%% Definition
Nph=500;  % photon number for 1p and 2p
K=9;
L=100;    %scan range in nanometers

%% Target coordinate position
vishift=-pi/8;

pos_nm(1,1)=L/2; pos_nm(1,2)=0;
pos_nm(1,3)=-vishift;
pos_nm(2,1)=L/2; pos_nm(2,2)=0;
pos_nm(2,3)=vishift;
pos_nm(3,1)=(L/2); pos_nm(3,2)=0;
pos_nm(3,3)=0;

faiorb=2*pi/3;
pos_nm(4,1)=(L/2)*cos(faiorb); pos_nm(4,2)=(L/2)*sin(faiorb);
pos_nm(4,3)=faiorb-vishift;
pos_nm(5,1)=(L/2)*cos(faiorb); pos_nm(5,2)=(L/2)*sin(faiorb);
pos_nm(5,3)=faiorb+vishift;
pos_nm(6,1)=(L/2)*cos(faiorb); pos_nm(6,2)=(L/2)*sin(faiorb);
pos_nm(6,3)=faiorb;

faiorb=4*pi/3;
pos_nm(7,1)=(L/2)*cos(faiorb); pos_nm(7,2)=(L/2)*sin(faiorb);
pos_nm(7,3)=faiorb-vishift;
pos_nm(8,1)=(L/2)*cos(faiorb); pos_nm(8,2)=(L/2)*sin(faiorb);
pos_nm(8,3)=faiorb+vishift;
pos_nm(9,1)=(L/2)*cos(faiorb); pos_nm(9,2)=(L/2)*sin(faiorb);
pos_nm(9,3)=faiorb;


%% Simulation
wavelength=640;
fwhm=wavelength/2/1.4;  
        
                     
       
Dntxyz = @(x,y,faiz) (1-cos(2*atan2(y,x)-2*faiz)).*4*exp(1)*log(2)*(x^2+y^2)/fwhm^2 .*exp(-4*log(2)*(x^2+y^2)/fwhm^2);
        
% zscale denote the scale in the axial direction    
zsca=wavelength/(2*pi*1.1);
I1 = @(x,y,faiz) Dntxyz(x-pos_nm(1,1),y-pos_nm(1,2),faiz/zsca+pos_nm(1,3));  
I2 = @(x,y,faiz) Dntxyz(x-pos_nm(2,1),y-pos_nm(2,2),faiz/zsca+pos_nm(2,3)); 
I3 = @(x,y,faiz) Dntxyz(x-pos_nm(3,1),y-pos_nm(3,2),faiz/zsca+pos_nm(3,3));
I4 = @(x,y,faiz) Dntxyz(x-pos_nm(4,1),y-pos_nm(4,2),faiz/zsca+pos_nm(4,3));  
I5 = @(x,y,faiz) Dntxyz(x-pos_nm(5,1),y-pos_nm(5,2),faiz/zsca+pos_nm(5,3)); 
I6 = @(x,y,faiz) Dntxyz(x-pos_nm(6,1),y-pos_nm(6,2),faiz/zsca+pos_nm(6,3));
I7 = @(x,y,faiz) Dntxyz(x-pos_nm(7,1),y-pos_nm(7,2),faiz/zsca+pos_nm(7,3));  
I8 = @(x,y,faiz) Dntxyz(x-pos_nm(8,1),y-pos_nm(8,2),faiz/zsca+pos_nm(8,3)); 
I9 = @(x,y,faiz) Dntxyz(x-pos_nm(9,1),y-pos_nm(9,2),faiz/zsca+pos_nm(9,3));

I_sum=@(x,y,faiz) I1(x,y,faiz)+I2(x,y,faiz)+I3(x,y,faiz)+I4(x,y,faiz)+I5(x,y,faiz)+I6(x,y,faiz)+I7(x,y,faiz)+I8(x,y,faiz)+I9(x,y,faiz);
I_sum0= I1(0,0,0)+I2(0,0,0)+I3(0,0,0)+I4(0,0,0)+I5(0,0,0)+I6(0,0,0)+I7(0,0,0)+I8(0,0,0)+I9(0,0,0);
    
%% SBR Definition
SBR=@(x,y,faiz,sbr0) (I_sum(x,y,faiz)/I_sum0)*sbr0;

%% MLE Definition
    P1 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I1(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P2 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I2(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P3 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I3(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P4 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I4(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P5 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I5(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P6 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I6(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P7 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I7(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P8 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I8(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    P9 = @(x,y,faiz,sbr0) (SBR(x,y,faiz,sbr0)/(SBR(x,y,faiz,sbr0)+1))*I9(x,y,faiz)/I_sum(x,y,faiz)+1/(K*(SBR(x,y,faiz,sbr0)+1));
    

    
    %% Simulation
    
%% define three molecules with position (0,0,0)nm, (15,15,15)nm and
% -(15,15,15)nm
ro(1,1)=0; r0(2,1)=0; r0(3,1)=0;
rdis=15; 
r0(1,2)=rdis;r0(2,2)=rdis; r0(3,2)=rdis;
r0(1,3)=-rdis;r0(2,3)=-rdis; r0(3,3)=-rdis;

times=1000; % number of simulated measurements
res=zeros(times*length(ro(1,:)),4);  
opt=optimset('Display','off');
res_avg_r0=[0,0,0,0];
prec_xy_r0=[0,0,0,0];

 %estimate for each molecules  
 for m=1:3  %length(ro(1,:))
    tic
     xx=r0(1,m);
     yy=r0(2,m);
     zz=r0(3,m);

    %change parfor to for here
    for jj=(1+(m-1)*times):(m*times)
        
        p1=P1(xx,yy,zz,5); p2=P2(xx,yy,zz,5); p3=P3(xx,yy,zz,5); 
        p4=P4(xx,yy,zz,5); p5=P5(xx,yy,zz,5); p6=P6(xx,yy,zz,5); 
        p7=P7(xx,yy,zz,5); p8=P8(xx,yy,zz,5); p9=P9(xx,yy,zz,5); 
        
        % generate random measurments based on the probalility and noise 
        R = mnrnd(Nph,[p1,p2,p3,p4,p5,p6,p7,p8,p9]);           
        %R=nexp(jj,:);  % nexp is the simulated measurements using true psf
        n1=R(1); n2=R(2); n3=R(3); 
        n4=R(4); n5=R(5); n6=R(6); 
        n7=R(7); n8=R(8); n9=R(9); 
        
        %four parameters are fitted, (x,y,z, sbr0)
        funcm=@(x)(-n1*log(P1(x(1),x(2),x(3),x(4)))+P1(x(1),x(2),x(3),x(4))-n2*log(P2(x(1),x(2),x(3),x(4)))+P2(x(1),x(2),x(3),x(4))-n3*log(P3(x(1),x(2),x(3),x(4)))+P3(x(1),x(2),x(3),x(4))...
            -n4*log(P4(x(1),x(2),x(3),x(4)))+P4(x(1),x(2),x(3),x(4))-n5*log(P5(x(1),x(2),x(3),x(4)))+P5(x(1),x(2),x(3),x(4))-n6*log(P6(x(1),x(2),x(3),x(4)))+P6(x(1),x(2),x(3),x(4))...
            -n7*log(P7(x(1),x(2),x(3),x(4)))+P7(x(1),x(2),x(3),x(4))-n8*log(P8(x(1),x(2),x(3),x(4)))+P8(x(1),x(2),x(3),x(4))-n9*log(P9(x(1),x(2),x(3),x(4)))+P9(x(1),x(2),x(3),x(4)) );
       
       % solve the molecule's position using fmincon 
       vlb=[-L;-L;-L; 1];
       vub=[L;L;L; 10];
       r_LMS1=[0,0,0, 3];
       [solution,fval,exitflag,output]=fmincon(@(x)funcm(x),r_LMS1,[],[],[],[],vlb,vub);
       disp(solution)

       res(jj,:)=solution;

    end
    res_avg_r0(m,:) = mean(res((1+(m-1)*times):(m*times),:)); % average
    prec_xy_r0(m,:) = std(res((1+(m-1)*times):(m*times),:));  % std (precision)
end
    
figure
for m=1:3
scatter3(res((1+(m-1)*times):(m*times),1),res((1+(m-1)*times):(m*times),2),res((1+(m-1)*times):(m*times),3),'filled')
hold on
end   
    
    
    