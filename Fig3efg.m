clear all; close all; clc

% DESCRIPTION - reproducing parts of Figure 3 e,f,g
% condenser NA is 0.9
% used functions:   off_axis_reconstruction_general 
%                   apodiz_SG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            LOADING DATA           %%%%%%%%%%%%%%%%%%%%%%

directory='...';  % path to the MAIN folder
directory_images=[directory '\data']; % folder with camera records

% direct image
Iint=single(imread([directory_images '\2024-11-06-1045_5-CapObj_0009.png']));

% interference pattern with yeast
I=apodization_for_propag(single(imread([directory_images '\2024-11-06-1046_1-CapObj_0002.png'])));      

% interference pattern without sample
Iref=apodization_for_propag(single(imread([directory_images '\2024-11-05-1219_2-CapObj_0000.png'])));

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

% position of the diffraction order in the Fourier spectrum
row_R=1392; col_R=2138;              % center of the usefull diffraction order in interf pattern with sample
row_R_ref=row_R; col_R_ref=col_R-1;  % center of the usefull diffraction order in ref interf. pattern

% Processing the interference pattern with sample
[U,~,~]=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R).*X+fy(row_R).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12);

% Processing the reference interference pattern
[Urefpom,~,~]=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R_ref).*X+fy(row_R_ref).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12);

% Corrected complex amplitude (corrected residual tilt and plane)
Uref=Urefpom./exp(1i*2*pi/lambda*((0.00004.*X)+(-0.00003*Y))).*exp(1i*0.7*pi);

% coordinates for plot
rows_1=(810:1210);
cols_1=(1550:2100);

unw_phase=angle(U(rows_1,cols_1)./Uref(rows_1,cols_1));
unw_phase=unw_phase-min(unw_phase(:)); %unwrapped phase

% Figure 3e
figure(1); 
imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,sqrt(Iint(rows_1,cols_1))); 
axis image; colormap gray; colorbar;
set(gca,'fontsize',15)
title('sqrt(I)')

% Figure 3f
figure(2); 
imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,abs(U(rows_1,cols_1))*mean2(sqrt(Iint(rows_1,cols_1)))/mean2(abs(U(rows_1,cols_1)))); 
axis image; colormap gray; colorbar;
set(gca,'fontsize',15); title('amplitude')

% Figure 3g
figure(3); 
imagesc(x(cols_1)/beta*1000000,y(rows_1)/beta*1000000,unw_phase); 
axis image; colormap parula; colorbar;
set(gca,'fontsize',15); title('phase')