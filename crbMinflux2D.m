
function [sigma_crb, sigma_crbx, sigma_crby]= crbMinflux2D(beadpsf, Nph, SBR)
%beadpsf: 3D excitation intensity profiles for multiple target coordinates 
%Nph: the photons detected 
%SBR: signal to background ratio in detection.

  K=size(beadpsf,3);
  step_nm=1;   %nm for x and y and z
  sizexy=size(beadpsf,1)-1;
  s=50;  %rough estimate precision

 
  %calculate the crlb
  d=2;   % d=2 for 2d
  dx=step_nm; dy=step_nm; 
  sigma_crb=zeros(sizexy,sizexy);  sigma_crbx=sigma_crb;  sigma_crby=sigma_crb; 
  Fr_aux=zeros(d,d,sizexy,sizexy,K);
  psf_sim=(beadpsf(:,:,1:K));
  sbr_rel=sum(psf_sim,3)./sum(psf_sim(sizexy/2+1,sizexy/2+1,:),3);
  sbr=sbr_rel.*SBR;
  norm_psfsim=sum(psf_sim,3);


  for i=1:K
    p(:,:,i)=(sbr./(sbr+1)).*psf_sim(:,:,i)./norm_psfsim+(1./(sbr+1)).*(1/K);
    [dpdx(:,:,i),dpdy(:,:,i)]=gradient(p(:,:,i),dx,dy);
  end

  prior='rough loc'
  for i=1:sizexy
        disp(i)
        for j=1:sizexy
              for k=1:K

                A=[dpdx(i,j,k)^2, dpdx(i,j,k)*dpdy(i,j,k);  
                    dpdx(i,j,k)*dpdy(i,j,k), dpdy(i,j,k)*dpdy(i,j,k)];
                    Fr_aux(:,:,i,j,k)=(1./p(i,j,k)).*A;
                    
               end
                if prior=='rough loc'
                    Fr(:,:,i,j)=Nph.*sum(Fr_aux(:,:,i,j,:),5)+(1/s^2).* diag([1,1]);
                else 
                    Fr(:,:,i,j)=Nph.*sum(Fr_aux(:,:,i,j,:),5);
                end
               sum_crb(:,:,i,j)=inv(Fr(:,:,i,j));
               sigma_crb(i,j)=sqrt((1/d).*trace(sum_crb(:,:,i,j)));
               sigma_crbx(i,j) = sqrt((sum_crb(1,1,i,j)));
               sigma_crby(i,j) = sqrt((sum_crb(2,2,i,j)));

        end
   end


end

