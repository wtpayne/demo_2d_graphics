function gabor = crsCreateGaborPrimitive(output_size,location,gauss_angle,deviations,sine_angle,frequency,phase);
% CRSCREATEGABORPRIMITIVE returns a real valued Gabor function in the range [-1,1].
% -----------------------------------------------------------------------------
% crsCreateGaborPrimitive
% =======================
% This function is used to generate a single Gabor function with values
% ranging from -1 to +1, "background" value 0. (The "background" value and
% the mean value are not necessarily the same.)
%
% A Gabor function (Gabor wavelet) is commonly used as a model of the receptive
% field of simple cells, and stimuli constructed using Gabor wavelets have the
% neat property of finite extent in both the spatial and Fourier domains -
% potentially limiting spurious excitations.
%
% A Gabor function - however it is specified, is essentially a sine function
% modulated (multiplied) by a Gaussian envelope to limit its spatial extent.
%
% This is an auxiliary function, and it is intended that the returned 2D matrix
% should be scaled, translated to a given location on the screen, added or
% otherwise combined with other Gabors, then scaled to a range of pixel levels
% to provide the desired stimulus.
%
% Other CRS matlab functions will use this function internally, and expose a
% more user friendly interface: In particular one that allows the location - in
% cartesian coordinates - to be specified explicitly.
%
% If you do wish to use this function directly, it may be an idea to have a
% look at the matlab source code: "edit GaborFun" to get an idea of how the
% various parameters are used to generate the function.
%
% Usage
% =====
% gabor = ...
%   crsCreateGaborPrimitive( ...
%               output_size, ...
%               location,    ...
%               gauss_angle, ...
%               deviations,  ...
%               sine_angle,  ...
%               frequency,   ...
%               phase);
%
% Parameters
% ==========
% output_size - Measured in pixels, this parameter should be a length 2 vector,
%               containing the following values: [height,width]. height and
%               width should be positive integer values.
%
%               If you want the output to fill the entire screen, you can use
%               crsGetScreenSizePixels to get the correct value.
%
% location    - Measured in pixels, this parameter should be a length 2 vector,
%               containing the following values: [x,y]. x and y should both be
%               real values.
%
%               The location is measured from the centre of the output matrix
%               to the centre of the gaussian patch, and can be specified in
%               subpixel increments.
%
% gauss_angle - Measured in degrees, this parameter should be scalar, real and
%               lie between 0 and 360.
%
%               The gaussian envelope can be rotated to an arbitrary angle,
%               specified (in degrees) with this parameter, which should be a
%               scalar real value, lying between 0 and 360. At 0 degrees, the
%               "major" axis is horizontal.
%
% deviations  - Measured in samples, this parameter should be a length 2 vector
%               containing the following values: [major,minor]. major and minor
%               should be positive real values.
%
%               The 2D gaussian envelope is actually constructed by multiplying
%               together two orthgonal 1D gaussian distributions. By changing
%               the relative values of the two standard deviations, different
%               aspect ratios are possible - allowing gabor functions with
%               elliptical spatial extent as well as varying size
%
% sine_angle  - Measured in degrees, this parameter should be scalar, real and
%               lie between 0 and 360.
%
% frequency   - Measured in cycles per sample, this parameter should be scalar,
%               real, and positive. (It will also normally lie between 0 and 1)
%
% phase       - Measured in degrees, this parameter should be scalar, real and
%               positive. The phase is measured relative to the centre of the
%               Gaussian envelope.
%
% Return values
% =============
% gabor       - A matrix of size output_size, containing the generated gabor
%               function.
%
% Reference page in Help browser
% <a href="matlab:web(['jar:file:',which('crsHelp\help.jar'),'!/crs\tools\stimulusmatrices\crsCreateGaborPrimitive.html'],'-helpbrowser');">crsCreateGaborPrimitive HTML help.</a>                                                                                                                         
%
% -----------------------------------------------------------------------------

% TYPE CHECKING.
% ==============

[height,width] = size(output_size);
if (height~=2 && width~=2) || any(any(output_size<0))
  error('output_size must be an array of vectors of length 2 holding positive integers.');
end
if height == 2
  output_size = output_size';
end

[height,width] = size(location);
if (height~=2 && width~=2)
  error('location must be an array of vectors of length 2 holding real values.');
end
if height == 2
  location = location';
end

if min(size(gauss_angle)) ~= 1 || any(any(gauss_angle<0)) || any(any(360<gauss_angle))
  error('gauss_angle must be an array of scalar real values between 0 and 360.');
end

[height,width] = size(deviations);
if (height~=2 && width~=2)  || any(any(deviations<0))
  error('deviations must be an array of vectors of length 2 holding positive real values.');
end
if height == 2
  deviations = deviations';
end

if min(size(sine_angle)) ~= 1 || any(any(sine_angle<0)) || any(any(360<sine_angle))
  error('sine_angle must be an array of scalar real values between 0 and 360.')
end

if min(size(frequency)) ~= 1 || any(any(frequency<0))
  error('frequency must be an array of positive scalar values.')
end

if min(size(phase)) ~= 1
  error('phase must be a scalar real value.')
end

if length(phase) ~= length(frequency) || (length(phase) ~= length(sine_angle)) || length(phase) ~= length(gauss_angle) || length(location) ~= length(deviations)
  error('Parameters must have consistent sizes.');
end

% Create the sampling grid used for both Gaussians.
% =================================================
% Allow at least 4 standard deviations on each side if we have a 15 bit display
% The value of the gaussian at the clip point will be insignificant.
max_spatial_extent = max(max(deviations)) * 8;

[xSampleGrid,ySampleGrid] = ...
ndgrid((max_spatial_extent-1)/2:-1:-((max_spatial_extent-1)/2), ...
       (max_spatial_extent-1)/2:-1:-((max_spatial_extent-1)/2));
clear max_spatial_extent;

% -----------------------------------------------------------------------------

gabor = zeros(output_size);
count_gabors = length(phase);

for i = 1:count_gabors

  % The gaussians are orthogonal (at 90 degrees to one another).
  % The means are kept at 0 (translations occur later).

  % First guassian.
  % ===============
  [xComp,yComp] = pol2cart(gauss_angle(i) * (pi/180),1);
  sampleGrid = (xSampleGrid * xComp) + (ySampleGrid * yComp);
  temp = exp((-sampleGrid.^2)/(2*(deviations(i,1)^2)));

  % Second guassian.
  % ================
  [xComp,yComp] = pol2cart((gauss_angle(i)+90) * (pi/180),1);
  sampleGrid = (xSampleGrid * xComp) + (ySampleGrid * yComp);
  temp = temp .* exp((-sampleGrid.^2)/(2*(deviations(i,2)^2)));

  % Sinusoid.
  % =========
  [xComp,yComp] = pol2cart(sine_angle(i) * (pi/180),1);
  sampleGrid    = (xSampleGrid * xComp) + (ySampleGrid * yComp);
  temp          = temp .* sin(phase(i)*(pi/180)+(frequency(i)*sampleGrid*2*pi));

  % Translate the gaussian to its final location
  % ============================================
  temp = crsTranslateMatrix(temp,location(i,:),output_size,0);

  % Accumulate Gabors
  % =================
  gabor = gabor + temp;

end

clear xComp yComp xSampleGrid ySampleGrid;
% -----------------------------------------------------------------------------
