clear all;  
close all;

% parameters 
% na: numerical aperture,  lambd: wavelength of the laser
% spixel: pixel size for excitation psf  
% num_x : size of the excitation psf
% LL 

%% basic parameters for calculating the excitation point spread function
na=1.33;      % numerical aperture of the optical system  
lambda=0.64;  % micrometer
%z=(0:0.1:1); % z's dimension should be 1*num_z1 
d=50; d1=[]; d2=[];
pix=0.06;%  micrometer  
mag=60;
spix=pix/mag;  % pixel size 
num_x=2*100+1; % number of pixels in the excitation point spread function 



%% define the coordinats in the fourier plane and sample plane
n1=1.52; n0=1.33;   %  refractive index for oil and sample
k0=2*pi/lambda; fpupil=[];
kmax=k0*na;
[kx_c xx_c]=czt_getcztcoord(num_x,300);
kx=(k0*na).*kx_c; 
xx=(num_x*spix/2).*xx_c;
[kxx1 kyy1]=meshgrid(kx);
[pxx1 pyy1]=meshgrid(xx);
NN=length(kx); MM=length(xx);
R=(k0*na) ;   L=(num_x*spix/2);
CTF1=circ(kxx1,kyy1,2*k0*na); 
%kz in the immersion oil
kz_c1=CTF1.*sqrt((k0*n1)^2-kxx1.^2-kyy1.^2);
%kz in the medium
kz_c2=CTF1.*sqrt((k0*n0)^2-kxx1.^2-kyy1.^2);
% define the apodization function for calculating the excitation point
% spread function
do_apod=1;
if do_apod==1
apod=sqrt(kz_c1/(n1.*k0));
else 
    apod=ones(nn,nn);
end


%% define the target coordinates in viminflux
vishift=pi/4;  LL=0.1;  %
K=10;
% (0,0) position is added to avoid ambiguity points, but it is not
% necessary in real acquisitions
pos_nm=zeros(K,3);
thetaorbit=zeros(K,1);
pos_nm(1,1)=0; pos_nm(1,2)=0; 
%fai shift
pos_nm(1,3)=0;

% TC1 and three phase shifted pattern
pos_nm(2,1)=LL/2; pos_nm(2,2)=0;
pos_nm(2,3)=-vishift;
pos_nm(3,1)=LL/2; pos_nm(3,2)=0;
pos_nm(3,3)=vishift;
pos_nm(4,1)=(LL/2); pos_nm(4,2)=0;
pos_nm(4,3)=0;

faiorb=2*pi/3;     
%faiorb=pi/2;
pos_nm(5,1)=(LL/2)*cos(faiorb); pos_nm(5,2)=(LL/2)*sin(faiorb);
pos_nm(5,3)=-mod(2*faiorb,2*pi)-vishift;
pos_nm(6,1)=(LL/2)*cos(faiorb); pos_nm(6,2)=(LL/2)*sin(faiorb);
pos_nm(6,3)=-mod(2*faiorb,2*pi)+vishift;
pos_nm(7,1)=(LL/2)*cos(faiorb); pos_nm(7,2)=(LL/2)*sin(faiorb);
pos_nm(7,3)=-mod(2*faiorb,2*pi);

faiorb=4*pi/3;
pos_nm(8,1)=(LL/2)*cos(faiorb); pos_nm(8,2)=(LL/2)*sin(faiorb);
pos_nm(8,3)=-mod(2*faiorb,2*pi)-vishift;
pos_nm(9,1)=(LL/2)*cos(faiorb); pos_nm(9,2)=(LL/2)*sin(faiorb);
pos_nm(9,3)=-mod(2*faiorb,2*pi)+vishift;
pos_nm(10,1)=(LL/2)*cos(faiorb); pos_nm(10,2)=(LL/2)*sin(faiorb);
pos_nm(10,3)=-mod(2*faiorb,2*pi);



%% calculate the 3D psfs using czt transform
x0=0; y0=0; 
fai=atan2(kxx1,kyy1);

for k=1:K
  fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0+pos_nm(k,1),y0+pos_nm(k,2),0,0);
  fmask1=exp(1i.*fai);         
  fkxky=fpupil.*fmask1.*apod;
  im=czt_getcztfft(fkxky,R,L,NN,MM);
  fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0+pos_nm(k,1),y0+pos_nm(k,2),0,0);
  fmask1=exp(-1i.*fai);
  fkxky=fpupil.*fmask1.*apod;
  im1=czt_getcztfft(fkxky,R,L,NN,MM);
  psfabc=abs(im+exp(1i*pos_nm(k,3)).*im1).^2;
  beadpsf(:,:,k)=psfabc;  % save beadpsf
 % figure(10)      % show the images when debugging
  %imagesc(psfabc);

end



% normalize the psf by the psf in the middle axial plane
sumpsf=sum(sum(beadpsf,1),2);
Imidpsf=sumpsf(1,1,1);   
beadpsf=1000.*beadpsf./Imidpsf;
    

% calculate the localization precision for 3D viMinflux
Nph=500; SBR=5;    
[sigma_crb, sigma_crbx, sigma_crby]= crbMinflux2D(beadpsf, Nph, SBR);

%define 2d as the mean of sigma x and sigma y
sigma_xy= (sigma_crbx+sigma_crby)/2;
roix=floor(1000*L/2); midx=101;

%plot the 2d precisions inside the scan range (-L/2: :L/2)
figure
imagesc(sigma_xy(midx-roix:midx+roix,midx-roix:midx+roix).*sqrt(Nph), [0 120]);
figure
plot(sigma_xy(midx-roix:midx+roix,midx).*sqrt(Nph))


