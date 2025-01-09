wl=0.7;  %wavelength of light in um
w=2*pi/6.8;   %spatial modulation sin(w*x)
n=1;  %refractive index of air
A=0.6;    %amplitude of modulation
rad=20; %radius of beam or sample illuminated 
if_aper=1; %if aperture-circular aperture for circular spot
per=0.2; %spatial resolution 
x_i=1.5*10^4;  %image plane distance in um

r_i=floor(rad/per);
x=per*linspace(-r_i,r_i,2*r_i+1);
y=per*linspace(-r_i,r_i,2*r_i+1);

k=2*pi/wl;
phase1=repmat(2*k*(n)*A.*sin(w.*x),[2*r_i+1,1]);
% phase1=repmat(-k.*(sqrt(x.^2+f^2)-f),[2*r_i+1,1]);

figure(1);
imagesc(phase1);
colorbar;


phi_0=0;%angle for incident plane
phi_0=phi_0*pi/180;

i0=35;%angle of incidence
i0=i0*pi/180;

phase2=k.*(x.*cos(phi_0)+y'*sin(phi_0)).*sin(i0);
% phase2=0;

figure(2)
image(x,y,phase2)
colorbar
% phase2(101,101)

phase=phase1+phase2;
aper=sqrt(x.^2+y'.^2);
aper(aper<rad)=1;
aper(aper>=rad)=0;
figure(3);
imagesc(aper);
colorbar;


if(if_aper==1)
    phase=phase.*aper;
end

figure(4);
imagesc(phase.*aper);
colorbar

E0=exp(1i.*phase);

% figure(5);
% imagesc(x,y,abs(E0).^2); colormap jet; axis xy;
% colorbar
% xlabel('x(um)');
% ylabel('y(um');
% title('object Plane');

res=500;   %no of points to calculate E field.More res means more clarity of image but more time.
cenz=x_i*cot(i0)*sec(phi_0);ceny=x_i*tan(phi_0);
y_im=linspace(-1.5*10^4+ceny,1.5*10^4+ceny,res);     %tuning y abd z image coordinates.ceny and cenz are pos of the reflected wave.
z_im=linspace(0,3.5*10^4+(cenz-1.5*10^4),res);
E_im=zeros(res,res);

for i=1:res
    for j=1:res
        r12=sqrt( (x-x_i).^2+(y'-y_im(i)).^2);
        l=sqrt(z_im(j)^2+r12.^2);
        E_im(i,j)=sum( (1/wl).*E0.*((1./(k.*l)-1i)).*z_im(j).*exp(1i.*k.*l)./l.^2,'all')*(per/2)^2;
    end
end
   
figure(6);
imagesc(y_im,z_im,abs(E_im').^2); axis xy;
colorbar
xlabel('y(um)');
ylabel('z(um');
title('image Plane');

