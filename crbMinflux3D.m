

function [sigma_crb, sigma_crbx, sigma_crby, sigma_crbz]= crbMinflux3D(beadpsf, Nph, SBR)
%beadpsf: 3D excitation intensity profiles for multiple target coordinates 
%Nph: the photons detected 
%SBR: signal to background ratio in detection.


  K=size(beadpsf,4);
  step_nm=1;   %nm for x and y and z
  sizexy=size(beadpsf,1)-1;
  sizez=size(beadpsf,3)-1;
  s=50;  %rough estimate precision

 
  %calculate the crlb
  d=3;
  dx=step_nm; dy=step_nm; dz=step_nm;
  sigma_crb=zeros(sizexy,sizexy, sizez);  sigma_crbx=sigma_crb;
  sigma_crby=sigma_crb; sigma_crbz=sigma_crb;
  Fr_aux=zeros(d,d,sizexy,sizexy,sizez,K);
  psf_sim=(beadpsf(:,:,:,1:K));
  sbr_rel=sum(psf_sim,4)./sum(psf_sim(sizexy/2+1,sizexy/2+1,sizez/2+1,:),4);
  sbr=sbr_rel.*SBR;
  norm_psfsim=sum(psf_sim,4);


  for i=1:K
    p(:,:,:,i)=(sbr./(sbr+1)).*psf_sim(:,:,:,i)./norm_psfsim+(1./(sbr+1)).*(1/K);
    [dpdx(:,:,:,i),dpdy(:,:,:,i), dpdz(:,:,:,i)]=gradient(p(:,:,:,i),dx,dy,dz);
  end

  prior='rough loc'
    for i=1:sizexy
        disp(i)
        for j=1:sizexy
           for m=1:sizez
              for k=1:K

                A=[dpdx(i,j,m,k)^2, dpdx(i,j,m,k)*dpdy(i,j,m,k),  dpdx(i,j,m,k)*dpdz(i,j,m,k);  
                    dpdx(i,j,m,k)*dpdy(i,j,m,k), dpdy(i,j,m,k)*dpdy(i,j,m,k), dpdy(i,j,m,k)*dpdz(i,j,m,k);
                    dpdx(i,j,m,k)*dpdz(i,j,m,k),  dpdy(i,j,m,k)*dpdz(i,j,m,k), dpdz(i,j,m,k).^2];
                Fr_aux(:,:,i,j,m,k)=(1./p(i,j,m,k)).*A;
               end
                if prior=='rough loc'
                    Fr(:,:,i,j,m)=Nph.*sum(Fr_aux(:,:,i,j,m,:),6)+(1/s^2).* diag([1,1,1]);
                else 
                    Fr(:,:,i,j,m)=Nph.*sum(Fr_aux(:,:,i,j,m,:),6);
                end
               sum_crb(:,:,i,j,m)=inv(Fr(:,:,i,j,m));
               sigma_crb(i,j,m)=sqrt((1/d).*trace(sum_crb(:,:,i,j,m)));
               sigma_crbx(i,j,m) = sqrt((sum_crb(1,1,i,j,m)));
               sigma_crby(i,j,m) = sqrt((sum_crb(2,2,i,j,m)));
               sigma_crbz(i,j,m) = sqrt((sum_crb(3,3,i,j,m)));

           end
        end
    end

end

