function IMAGE = crsCreateFourierNoisePrimitive(IMG_SIZE,FUN);
% CRSCREATEFOURIERNOISEPRIMITIVE returns noise with a specified power spectrum.
% -----------------------------------------------------------------------------
% crsCreateFourierNoisePrimitive
% ==============================
% The crsCreateFourierNoisePrimitive function allows you to generate 
% anisotropic monochrome "noise" images with a user specified power spectrum.
%
% The images that are returned are two-dimensional matrices (i.e. monochrome
% images), with values in the range 0-1.
%
% If you want to view the images with MATLAB's image function, you will need to
% scale the image to the range 1-65, (depending on the length of the current
% colormap).
%
% The range 0-1 is suitable for combining stimuli through multiplication. If
% you want to combine stimuli with addition or through some other means, you
% will probably need to re-scale the resulting stimulus image before drawing it
% to video memory.
%
% Note: Stimuli matrices can be copied to VSG video memory using the
% crsDrawMatrix* functions. Please refer to the documentation for each of these
% functions to find out what their preferred range of input values is,
% (together with the appropriate video mode).
%
% You need to specify a power spectrum for this function. A flat power
% spectrum: ( ones(1,10), for example. ) gives white noise, whereas a 1/f
% power spectrum gives "pink" noise. A more extreme shift of power to
% lower frequencies is illustrated in the spectrum below:
%
% Usage
% =====
% FUN = fliplr(([1:200]/200).^20); % DC at start of vector.
% IMG_SIZE = 600;
% IMG = crsCreateFourierNoisePrimitive(IMG_SIZE,FUN);
% imagesc(IMG);
% crsDrawMatrix(IMG);
%
% If you want to create noise with isotropic characteristics, it should be
% possible to modify the source description in this file in an appropriate
% manner. Creating colour noise may be more difficult, although creating
% a three dimensional noise field that includes a temporal dimension (a
% noise movie) may actually be quite easy (although computationally demanding).
%
% Reference page in Help browser
% <a href="matlab:web(['jar:file:',which('crsHelp\help.jar'),'!/crs\tools\stimulusmatrices\crsCreateFourierNoisePrimitive.html'],'-helpbrowser');">crsCreateFourierNoisePrimitive HTML help.</a>                                                                                                           
%
% -----------------------------------------------------------------------------

FUN_SIZE = numel(FUN);

% Create a polar coordinate space., with
% the origin in the centre
[XGRID,YGRID] = meshgrid([0:IMG_SIZE]-(IMG_SIZE/2),[0:IMG_SIZE]-(IMG_SIZE/2));
[TH,RAD]      = cart2pol(XGRID,YGRID);

% Linearly map the radial component onto
% the range required to index into the
% supplied function.
RAD = RAD - min(min(RAD));
RAD = RAD ./ max(max(RAD)) .* (FUN_SIZE-1);
RAD = RAD + 1;

% Use the radial component to index into the
% fourier-domain envelope, then shift the
% whole thing so that the first element of
% the envelope corresponds to the DC component.
RAD_FUN     = interp1(FUN,RAD);
RAD_FUN     = ifftshift(RAD_FUN);

% Create some random coefficients in the fourier domain.
PHASE_RAND  = rand(size(RAD_FUN)) .* 2 .* pi;
AMPL_RAND   = rand(size(RAD_FUN)) .* RAD_FUN; % modulated by the radial envelope

% Convert these to a complex representation,
% and do the inverse fourier transform.
[REAL,IMAG] = pol2cart(PHASE_RAND,AMPL_RAND);
CMPLX       = complex(REAL,IMAG);
IMAGE       = ifft2( CMPLX );

% We have a complex image as we did not bother ensuring
% that the random coefficients were symmetric. No bother,
% just arbitrarily discard one of the components. (The
% imaginary component in this case).
IMAGE       = real( IMAGE );

% Map the image to a sensible range of values (0-1).
IMAGE       = (IMAGE - min(min(IMAGE)));
IMAGE       = ( IMAGE ./ max(max(IMAGE)) );

