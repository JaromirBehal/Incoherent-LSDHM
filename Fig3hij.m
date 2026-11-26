clear all; close all; clc

% DESCRIPTION - reproducing parts of Figure 3 h,i,j
% condenser NA is 1.4
% used functions:   off_axis_reconstruction_general 
%                   apodiz_SG
%                   apodization_for_propag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            LOADING DATA           %%%%%%%%%%%%%%%%%%%%%%

directory='...';  % path to the MAIN folder
directory_images=[directory '\data']; % folder with camera records

% direct image
Iint=(single(imread([directory_images '\2024-11-08-1157_0-CapObj_0003.png'])));

% interference pattern with yeast
I=apodization_for_propag(single(imread([directory_images '\2024-11-08-1156_3-CapObj_0004.png']))); 

% interference pattern without sample
Iref=apodization_for_propag(single(imread([directory_images '\2024-11-08-0929_8-CapObj_0003.png'])));

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
    NA=1.4;                 % numerical aperture of Microscope Objective (MO)
    magMO=100;              % nominal magnification of the MO
    fTL_ideal=0.180;        % focal length of Tube lens the MO is designed for
    fTL=0.300;              % focal length of Tube lens used in the experiment
    beta=magMO*fTL/fTL_ideal*(150/200); % overal lateral magnification of sample-imaging system    

row_R=614; col_R=1043; % center of the usefull diffraction order 

% Processing the interference pattern with sample
[U,~,~]=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R).*X+fy(row_R).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12);

% Processing the reference interference pattern
[Uref,~,~]=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R).*X+fy(row_R).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12);

% coordinates for plot
rows_1=(380:780); 
cols_1=(1250:1900); 

unw_phase=angle(U(rows_1,cols_1)./Uref(rows_1,cols_1))+1.1;  % phase    

% SUPPLEMENTARY FILE FIGURE S3 - wrapped phase
figure(1); 
subplot(1,3,1); imagesc(x/beta*10^6,y/beta*10^6,angle(U)); axis image; colormap parula; clim([-pi pi]);
subplot(1,3,2); imagesc(x/beta*10^6,y/beta*10^6,angle(Uref)); axis image; colormap parula; clim([-pi pi]);
subplot(1,3,3); imagesc(x/beta*10^6,y/beta*10^6,angle(U./Uref)); axis image; colormap parula; clim([-pi pi]);

% Figure 3h
figure(2); 
imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,sqrt(Iint(rows_1,cols_1))); 
axis image; colormap gray; colorbar;
set(gca,'fontsize',15); title('sqrt(I)')

% Figure 3i
figure(3); 
imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,abs(U(rows_1,cols_1))*mean2(sqrt(Iint(rows_1,cols_1)))/mean2(abs(U(rows_1,cols_1)))); 
axis image; colormap gray; colorbar;
set(gca,'fontsize',15); title('amplitude')

% Figure 3j
figure(4); 
imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,unw_phase); 
axis image; colormap parula; colorbar; clim([0 4.1])
set(gca,'fontsize',15); title('parula')