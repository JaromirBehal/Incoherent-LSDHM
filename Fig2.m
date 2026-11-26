clear all; close all; clc

% DESCRIPTION - reproducing parts of Figure 2
% condenser NA is 0.4
% used functions:   apodization_for_propag
%                   off_axis_reconstruction_general 
%                   apodiz_SG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            LOADING DATA           %%%%%%%%%%%%%%%%%%%%%%

directory='...'; % path to the MAIN folder
directory_images=[directory '\data']; % folder with camera records

% direct image
I_USAF=single(imread([directory_images '\USAF_1951_HR_direct.png']));

% interference pattern with cheek cells
I=apodization_for_propag(single(fitsread([directory_images '\Cheek_4.FIT'])));

% interference pattern without sample
Iref=apodization_for_propag(single(fitsread([directory_images '\Diatoms_50_ref.FIT'])));

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
    NA=0.25;                % numerical aperture of Microscope Objective (MO)
    magMO=10;               % nominal magnification of the MO
    fTL_ideal=0.180;        % focal length of Tube lens the MO is designed for
    fTL=0.300;              % focal length of Tube lens used in the experiment
    beta=magMO*fTL/fTL_ideal*(150/200); % overal lateral magnification of sample-imaging system   

row_R=472;     col_R=2355;      % center of the usefull diffraction order in interf pattern with sample
row_R_ref=469; col_R_ref=2359;  % center of the usefull diffraction order in ref interf. pattern

% Processing the interference pattern with sample
[U,~,~]=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R).*X+fy(row_R).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12);

% Processing the reference interference pattern
[Uref,~,~]=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R_ref).*X+fy(row_R_ref).*Y)),fx,fy,2.0*(NA/lambda)/beta,0.5,12);

% Corrected complex amplitude
Urec=(U./Uref.*exp(-1i*k*(0.00002*X))).*exp(-1i*0.92);

% coordinates of cross section in direct image of USAF 1951
x_USAF = [x(1479) x(1570)]/beta*1000000;
y_USAF = [y(1136) y(993)]/beta*1000000;


% Figure 2a
figure(1); imagesc(x/beta*1000000,y/beta*1000000,I_USAF); axis image; colormap gray; colorbar; xlim([-25 60]); ylim([-35 55])
% Figure 2b
figure(2); improfile(x/beta*1000000,y/beta*1000000,I_USAF,x_USAF,y_USAF); grid on;
% Figure 2c
figure(3); imagesc(x/beta*1000000,y/beta*1000000,I); axis image; colormap gray; colorbar;
% Figure 2d
figure(4); imagesc(fx*beta/1000,fy*beta/1000,log1p(abs(fftshift(fft2(I))))); axis image; colormap gray; clim([9 22]); colorbar 
% Figure 2e
figure(5); imagesc(x/beta*1000000,y/beta*1000000,abs(U)); axis image; colormap gray;
% Figure 2f
figure(6); imagesc(angle(Urec)); axis image; colormap parula; colorbar;