%%
% load different component Ex, Ey, Ez of field computed using
% Lumerical-FDTD
load('Ex_Airy_grating_11p_normalIncidence_tfsf_600_900_p745_w200.mat');
load('Ey_Airy_grating_11p_normalIncidence_tfsf_600_900_p745_w200.mat');
load('Ez_Airy_grating_11p_normalIncidence_tfsf_600_900_p745_w200.mat');

%%
% concotanate the different field component in one electric field E

E = zeros(size(Ex,1), size(Ex,2), 3, size(Ex,4));
E(:,:,1,:) = Ex;
E(:,:,2,:) = Ey;
E(:,:,3,:) = Ez;

%%
%%%%%%%%%%%%%%% dispersion relation for plasmons %%%%%%%%%%%%%%%%%
%wavelength of interest; defined by user
wl=0.800;
% epsilon air
epsd=1;
% epsilon Gold 'Johnson & Christy'
load('eps_Johnson_Christy.mat');
% calculate the dispersion relation for the plasmons
k0_interp=2*pi./LAM.*sqrt((epsd.*e_interp)./(epsd+e_interp));
% wavevector at discrete wavelength wl = 0.800µm
k0 = k0_interp(LAM==wl);

%%
%%%%%%%%%%%%% FT- Beam Propagation Method %%%%%%%%%%%%%%%%%%%%%%%%
% % padding of zeros to get better freq resolution
pad = 1000;
% numerical spatial resolution
res=0.01;
% initial plane for which we calculate forward and backward propagation
z0 = 9.1;     
% maximuam propagation distance
z1 = 19.1;  
% incremental step-size in the propagation direction 'z'
delta_z = linspace(0, 10, 201);
% wavlength ouver which Lumerical fields are simulated
wavelengths = 0.6:0.01:0.9;
pxl_wl = find(wavelengths==wl);

for p=1:3
    E0=squeeze(E(:,round(z0/res),p,pxl_wl));
    for loop=1:length(delta_z)
        E_prop(:,loop,p)=propagator(E0,k0,delta_z(loop),pad,res);
    end
end

%%
% intensity of the fields calculated using FT-BPM
INT_BPM=abs(squeeze(E_prop(pad+1:end-pad,:,1))).^2+abs(squeeze(E_prop(pad+1:end-pad,:,2))).^2+...
            abs(squeeze(E_prop(pad+1:end-pad,:,3))).^2;

INT_Lumerical=abs(squeeze(E(:,:,1,pxl_wl))).^2+abs(squeeze(E(:,:,2,pxl_wl))).^2+...
            abs(squeeze(E(:,:,3,pxl_wl))).^2;
        

%%
% close all
x_Lumerical = linspace(0,30,3001);
y_Lumerical = linspace(0,10,1001);
figure(1);set(gcf,'Position',[500 50 900 600]);
h = pcolor(x_Lumerical,y_Lumerical,log10(INT_Lumerical(501:end-500,910:1910)).');
set(h, 'EdgeColor', 'none');
caxis([-2 1]);
cb = colorbar();
t = get(cb, 'Limits');
set(cb, 'Ticks', linspace(t(1), t(2), 3));
% ylabel(cb, 'log_{10} Intensity [a.u.]');
axis image;
% xlabel('x [µm]');
% ylabel('propagation distance [µm]');

set(gca,'xtick', [], 'ytick', [], 'fontsize',28, 'fontweight','bold','linewidth',2);

% print('-dpng','-r600','semiAnalytical_time_avg_intensity_25p_6fs_800nm.png')

%%
% close all
x_BPM = linspace(0,30,3001);
y_BPM = linspace(0,10,201);
figure(2);set(gcf,'Position',[500 50 900 600]);
h = pcolor(x_BPM,y_BPM,log10(INT_BPM(501:end-500,:)).');
set(h, 'EdgeColor', 'none');
caxis([-2 1]);
cb = colorbar();
t = get(cb, 'Limits');
set(cb, 'Ticks', linspace(t(1), t(2), 3));
% ylabel(cb, 'log_{10} Intensity [a.u.]');
axis image;
% xlabel('x [µm]');
% ylabel('propagation distance [µm]');

set(gca,'xtick', [], 'ytick', [], 'fontsize',28, 'fontweight','bold','linewidth',2);

% print('-dpng','-r600','semiAnalytical_time_avg_intensity_25p_6fs_800nm.png')


%%
I_Lumerical = zeros(3,3001);
I_Lumerical(1,:) = INT_Lumerical(501:end-500,910)/max(INT_Lumerical(501:end-500,910));
I_Lumerical(2,:) = INT_Lumerical(501:end-500,1410)/max(INT_Lumerical(501:end-500,1410));
I_Lumerical(3,:) = INT_Lumerical(501:end-500,1910)/max(INT_Lumerical(501:end-500,1910));

I_BPM = zeros(3,3001);
I_BPM(1,:) = INT_BPM(501:end-500,1)/max(INT_BPM(501:end-500,1));
I_BPM(2,:) = INT_BPM(501:end-500,101)/max(INT_BPM(501:end-500,101));
I_BPM(3,:) = INT_BPM(501:end-500,201)/max(INT_BPM(501:end-500,201));

%%
close all;
% z = 0µm
for loop = 1:3
    figure(100+loop); set(gcf,'Position',[300 50 900 600]);
    plot(x_Lumerical,I_Lumerical(loop,:),'--r','linewidth',1.5);
    hold on;
    plot(x_BPM,I_BPM(loop,:),'--b','linewidth',1.5);
    
    legend('Lumerical', 'FT-BPM');
    box on;
    grid off;
    xlabel('x [µm]');
    ylabel('intensity (norm.)');
    set(gca,'fontsize',16, 'fontweight','bold','linewidth',1.5);
    
%     print('-dpng','-r600',['intensity_lumerical_BPM_' num2str(5*(loop-1)) 'µm.png'])

end

