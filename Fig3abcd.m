clear all; close all; clc

% DESCRIPTION - reproducing parts of Figure 3 a,b,c,d
% condenser NA is 0.4
% used functions:   off_axis_reconstruction_general 
%                   apodiz_SG
%                   apodization_for_propag 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            LOADING DATA           %%%%%%%%%%%%%%%%%%%%%%

directory='...';  % path to the MAIN folder
directory_images=[directory '\data']; % folder with camera records

% direct image
Iint=apodization_for_propag(single(imread([directory_images '\2024-11-01-1442_9-CapObj_0000.png'])));

% interference pattern with diatom
I=apodization_for_propag(single(imread([directory_images '\2024-11-01-1441_4-CapObj_0000.png'])));

% interference pattern without sample
Iref=apodization_for_propag(single(imread([directory_images '\2024-11-01-1439_1-CapObj_0000.png'])));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%         OPTICAL PARAMETERS        %%%%%%%%%%%%%%%%%%%%%%
    % columns               % rows
    M=size(I,2);            N=size(I,1);
    m=-M/2:(M/2-1);         n=-N/2:(N/2-1);
    dx=2.4*10^(-6);         dy=dx;  % pixel size of the camera
    x=m.*dx;                y=n.*dy;
    [X,Y]=meshgrid(x,y);
    lambda=532*10^(-9); % wavelength    
    k=2*pi/lambda; % wavenumber

    fx=m/(M*dx);            fy=n/(N*dy); % spatial frequencies of the camera
    [FX,FY]=meshgrid(fx,fy);
    NA=0.4;                 % numerical aperture of Microscope Objective (MO)
    magMO=20;               % nominal magnification of the MO
    fTL_ideal=0.180;        % focal length of Tube lens the MO is designed for
    fTL=0.300;              % focal length of Tube lens used in the experiment
    beta=magMO*fTL/fTL_ideal*(150/200); % overal lateral magnification of sample-imaging system      

row_R=696; col_R=2151;              % center of the usefull diffraction order 
row_R_ref=row_R; col_R_ref=col_R;   % center of the usefull diffraction order 

% Processing the interference pattern with sample 
[U,~,~]=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R).*X+fy(row_R).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12); 

% Processing the reference interference pattern
[Uref,~,~]=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R_ref).*X+fy(row_R_ref).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12);

% For results in Fig. S2ac in Supplemental Document uncomment the two lines below
% [U,~,~]=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R).*X+fy(row_R).*Y)),fx,fy,1.0*(NA/lambda)/beta,0.5,12); 
% [Uref,~,~]=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R_ref).*X+fy(row_R_ref).*Y)),fx,fy,1.0*(NA/lambda)/beta,0.5,12);

% coordinates for plot
rows_1=1100:1560;
cols_1=790:1300;

% coordinates for improfile
x1=[339-90-15 410+70+15]+3; 
y1=[350+90+15 270-70-15]; 

unw_phase=(-angle(U(rows_1,cols_1)./Uref(rows_1,cols_1).*exp(+1i*2.2)))+2.4; % basic 'unwrapping' and shifting global phase
[cx,cy,c] = improfile(unw_phase,x1,y1);

dist_pix=sqrt((x1(1)-x1(2))^2+(y1(1)-y1(2))^2);
dist_real=(((1:length(c))/length(c))*dist_pix)*dx/beta*1000000; % linear distance of the improfile

% Figure 3a,b 
figure(1); 
subplot(1,2,1); imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,sqrt(Iint(rows_1,cols_1))); 
axis image; colormap gray; colorbar; clim([0 12]);
set(gca,'fontsize',15);
title('sqrt of intensity')
%
subplot(1,2,2); imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,abs(U(rows_1,cols_1))*mean2(sqrt(Iint(rows_1,cols_1)))/mean2(abs(U(rows_1,cols_1)))); 
axis image; colormap gray; colorbar; clim([0 12]);
set(gca,'fontsize',15);
title('amplitude')

% Figure 3c 
figure(2)
imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,unw_phase); axis image; colormap gray; colorbar; %clim([-pi pi])
hold on
plot((x([x1(1),x1(2)]+cols_1(1)))/beta*1000000,y(([y1(1),y1(2)]+rows_1(1)))/beta*1000000,'Color','r','LineWidth',2)
set(gca,'fontsize',15)
title('unwrapped phase')
colormap parula

% Figure 3d
figure(3)
plot(dist_real,c,'r','LineWidth',2); axis tight; ylim([-0.2 3.6])
title('profile')
set(gca,'fontsize',15)