function [] = ode_based_denoising(sigma,imnoisy,imoriginal)
Ts=200;
dt=.001;
if nargin==3
    flag=1;
else
    flag=0;
end
u0=imnoisy;
[mm,nn,channel]=size(imnoisy);
u_reconstructed=zeros(mm,nn,channel);
for j=1:channel
    u=imnoisy(:,:,j);
    Psnr= zeros(Ts+1,1);
    if flag==1
    Psnr(1) = psnr(imoriginal(:,:,j),uint8(u));
    else
        Psnr(1)=0;
        uold=imgaussfilt(u,2*sigma);
    end
for i = 1:Ts
    GU=imgaussfilt(u,sigma);
    upx = dpx(GU); umx = dmx(GU);  upy = dpy(GU); umy = dmy(GU);  ux = m(upx, umx);  uy = m(upy, umy);
    gradient_GU= sqrt(ux.^2 + uy.^2);
    M=max(max(gradient_GU(:)));
    u = u + dt*(gradient_GU.*(u-imgaussfilt(u,3*sigma))-M.*(u-imgaussfilt(u,sigma)));
if flag==1
    Psnr(i+1) = psnr(imoriginal(:,:,j),uint8(u));
else
    Psnr(i+1) = psnr(uint8(uold),uint8(u));
end
    if Psnr(i+1)<Psnr(i)
        break;
    end
end
u_reconstructed(:,:,j)=u;
end
if flag==1
 ssim(uint8(u_reconstructed(:,:,1)),imoriginal(:,:,1))
 psnr(imoriginal,uint8(u_reconstructed))
end
 figure
subplot(1,2,1)
imshow(uint8(u0))
title('Noisy Image')
subplot(1,2,2)
imshow(uint8(u_reconstructed))
title('Reconstructed Image')
return


%%%%%%%%%%%%%%%%%%%%%%%%%% minmod function
function c = m(a,b)
c = (a.*b>0).*sign(a).*min(abs(a),abs(b));
return
%%%%%%%%%%%%%%%%%%%%%%%%%%  forward x direction derivative function
function du = dpx(u)
du = u([2:end end],:) - u;
return
%%%%%%%%%%%%%%%%%%%%%%%%%% backward x direction derivative function
function du = dmx(u)
du = u - u([1 1:end-1],:);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%  forward y direction derivative function
function du = dpy(u)
du = u(:,[2:end end]) - u;
return
%%%%%%%%%%%%%%%%%%%%%%%%%% backward y direction derivative function
function du = dmy(u)
du = u - u(:,[1 1:end-1]);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%
