function gaussian_mesh = crsCreateGaussPrimitive(output_size,size_units,locations,locations_units,deviations,deviations_units,correlations,angles,scaling_factors);
% CRSCREATEGAUSSPRIMITIVE returns a sum of real valued gaussian functions
% -----------------------------------------------------------------------------
% crsCreateGaussPrimitive
% =======================
%
% The bivariate Gaussian function, or normal distribution, has the significant
% property that the Fourier transform of a Gaussian distribution is another
% (different) Gaussian distribution. The Gaussian distribution therefore has
% limited extent in both spatial and spatial-frequency domains.
%
% The bivariate Gaussian primitive appears as an elliptical patch on the image
% plane, smoothly varying from its perimiter to its center, with a simple
% Gaussian profile across all its diameters. The following function is used to
% generate the primitive: ( refer to the source for a complete explanation )
%
% ( k * exp(                           ...
%            -(   ((x1.^2)/sd1.^2)     ...
%               + ((x2.^2)/sd2.^2)     ...
%               - (   (2.*rho.*x1.*x2) ...
%                    / (sd1.*sd2)       ))));
%
% ( The above was adapted from the reference for a bivariate normal
%   distribution on http://mathworld.wolfram.com/ )
%
% The Gaussian profile ranges from a minimum of 0 to a maximum of 1. If more
% than one Gaussian is specified, then the Gaussians are added together, and
% the maximum value of the matrix returned may be greater than 1.
%
% The locations/means (in pixels), deviations(in pixels) and correllation
% coefficients can all be specified.
%
% This is an auxiliary function, and it is intended that the returned matrix
% should be scaled to a proportional value, multiplied or otherwise combined
% with other primitives, then scaled to a range of pixel levels to provide the
% desired stimulus before being transferred to the VSG using one of the
% drawing functions crsDrawMatrix or similar.
%
% Other CRS matlab functions will use this function internally, and expose a
% more user friendly interface: In particular one that allows units other than
% samples (pixels) to be used.
%
% If you do wish to use this function directly, it may be an idea to have a
% look at the matlab source code: "edit CreateGaussPrimitivePixels" to get
% an idea of how the various parameters are used to generate the function.
%
% Usage
% =====
% sum_of_gaussians = ...
%   crsCreateGaussPrimitive(  ...
%           output_size,      ...
%           size_units,       ...
%           locations,        ...
%           locations_units,  ...
%           deviations,       ...
%           deviations_units, ...
%           correlations,     ...
%           angles,           ...
%           scaling_factors);
%
% Parameters
% ==========
% output_size - Measured in either pixels, mm or degrees subtended at the
%               eye, this parameter should be a length 2 vector, containing
%               the following values: [height,width]. height and width should
%               be positive real values.
%
%               If you want the output to fill the entire screen, you can use
%               crsGetScreenSizePixels or crsGetScreenSizeDegrees to get the
%               correct value.
%
% size_units  - The units the output matrix size is measured in.
%               This can be one of the following constants:
%
%               CRS.MMUNIT     - mm across the stimulus display surface
%                                (monitor linewidth must be calibrated)
%
%               CRS.DEGREEUNIT - degrees subtended at the eye
%                                (monitor linewidth must be calibrated and
%                                 viewing distance must be specified
%                                 (default is 1000mm))
%
%               CRS.PIXELUNIT  - stimulus display device pixels
%
%
% locations   - Measured in either pixels, mm or degrees subtended at the
%               eye, this parameter should be a length 2 vector, containing
%               the following values: [x,y]. x and y should be real values.
%
%               The location is measured from the centre of the output matrix
%               to the centre of the gaussian patch, and can be specified in
%               subpixel increments.
%
%               If you want to obtain a sum of gaussians, then a vector of
%               locations may be provided. (giving a 2-by-n or n-by-2 matrix).
%               Note, however, that the same number of locations must be given
%               as deviations, correllations, angles and scaling_factors.
%
%               In the ambiguous case of a 2-by-2 matrix, we assume that
%               location vectors are held in rows.
%
% locations_units - The units that the locations vectors are measured in.
%                   This can be one of the following constants:
%
%                   CRS.MMUNIT     - mm across the stimulus display surface
%                                    (monitor linewidth must be calibrated)
%
%                   CRS.DEGREEUNIT - degrees subtended at the eye
%                                    (monitor linewidth must be calibrated and
%                                     viewing distance must be specified
%                                     (default is 1000mm))
%
%                   CRS.PIXELUNIT  - stimulus display device pixels
%
% deviations  - Measured in either pixels, mm, or degrees subtended at the eye,
%               this parameter should be a length 2 vector containing the
%               following values: [major,minor]. major and minor should be
%               positive real values.
%
%               For a patch where the functional form is recorded with 8-bit
%               numbers, the minimum representable value is just over 2.5
%               deviations from the centre. Where the functional form is
%               recorded with 14-bit numbers, this increases to just over 3
%               deviations.
%
%               The drop off in intensity at a range of standard deviations
%               is approximated by the following table:
%
%                              deviation      percentage
%                             from centre    of peak value
%                             ============   =============
%                                  0.0           100%
%                                  0.2            96%
%                                  0.4            85%
%                                  0.6            70%
%                                  0.8            54%
%                                  1.0            38%
%                                  1.2            24%
%                                  1.4            14%
%                                  1.6             8%
%                                  1.8             4%
%                                  2.0             1%
%
%               Another useful rule of thumb concerns the amount of "energy"
%               in the stimulus within a given number of deviations of the
%               centre:
%
%                 68%   of the energy should be within 1 s.d. of the center.
%                 95%   of the energy should be within 2 s.d. of the center.
%                 99.7% of the energy should be within 3 s.d. of the center.
%                 99.9% of the energy should be within 4 s.d. of the center.
%
%               In the ambiguous case of a 2-by-2 matrix, we assume that
%               deviation vectors are held in rows.
%
% deviations_units - The units that the deviations vectors are measured in.
%                    This can be one of the following constants:
%
%                    CRS.MMUNIT     - mm across the stimulus display surface
%                                     (monitor linewidth must be calibrated)
%
%                    CRS.DEGREEUNIT - degrees subtended at the eye
%                                     (monitor linewidth must be calibrated and
%                                      viewing distance must be specified
%                                      (default is 1000mm))
%
%                    CRS.PIXELUNIT  - stimulus display device pixels
%
% correlations - Lying in the interval [0,1), this parameter shoud be a vector
%                containing the same number of values as there are locations,
%                deviations, angles and scaling factors.
%
%                This parameter specifies the cross_correlation_coefficients
%                for each direction (x1,x2) in the output matrix.
%
%                A value close to 1 will cause the gaussian blob to extend
%                at 45 degrees between x1 and x2, forming a narrow extended
%                patch as the two axes become more highly correlated. A value
%                close to 0 will cause the gaussian blob to contract to a
%                symmetric, uncorrelated patch.
%
%                If you have a correlation close to 1, then the angle of the
%                resulting line-segment-like patch will be determined by
%                the relative values of the two deviation parameters for
%                this patch.
%
% angles       - Lying in the interval [0,2pi], this parameter specifies the
%                angle to which to rotate the axes containing each gaussian
%                patch. Strictly speaking, the same effect could be achieved
%                by manipulating the correlation and deviations parameters,
%                but using this parameter may be more intuitive for some.
%
%                This parameter has units of radians (of patch rotation), and
%                similarly to the other functions, the same number of angles
%                must be supplied as locations, deviations, correlations and
%                scaling factors.
%
% scaling_factors - This parameter specifies a scaling factor for each gaussian.
%                   This might be useful, for example, in constructing a DOG
%                   (Difference Of Gaussians) functional form - simply by scaling
%                   one of the gaussians by a negative number.
%
%                   Again, the same number of scaling_factors must be supplied
%                   as for the other vector parameters.
%
% Return values
% =============
% sum_of_gaussians - A matrix of the requested size, containing the sum of the
%                    generated gaussian functions.
%
% Example
% =======
%
%  output_size      = [600,800];
%  size_units       = CRS.PIXELUNIT;
%  locations        = [0,0];
%  locations_units  = CRS.PIXELUNIT;
%  deviations       = [50,50]
%  deviations_units = CRS.PIXELUNIT;
%  correlations     = 0;
%  angles           = 0;
%  scaling_factors  = 1;
%
% sum_of_gaussians = ...
%   crsCreateGaussPrimitive(  ...
%           output_size,      ...
%           size_units,       ...
%           locations,        ...
%           locations_units,  ...
%           deviations,       ...
%           deviations_units, ...
%           correlations,     ...
%           angles,           ...
%           scaling_factors);
%
% Reference page in Help browser
% <a href="matlab:web(['jar:file:',which('crsHelp\help.jar'),'!/crs\tools\stimulusmatrices\crsCreateGaussPrimitive.html'],'-helpbrowser');">crsCreateGaussPrimitive HTML help.</a>                                                                                                                         
%
% -----------------------------------------------------------------------------

  % We keep our constants in one place
  global CRS;

  % If we are given a single number for a size, the matrix is square
  output_size = round(crsUnitToUnit(size_units,output_size,CRS.PIXELUNIT));
  if length(output_size) == 1
    output_size = [output_size,output_size];
  end

  % Create a grid from -output_size/2 to +output_size/2
  if ( ~all(size(output_size) == [1,2]) && ~all(size(output_size) == [2,1]) )
    error('output_size must be a length 2 vector');
  end
  x1_vector = [0:(output_size(1)-1)] - (output_size(1)/2);
  x2_vector = [0:(output_size(2)-1)] - (output_size(2)/2);

  % we negate y to compensate for matrices being upside down :-)
  [x1_grid,x2_grid] = meshgrid(x1_vector,-x2_vector);

  % Calculate how many gaussians we have
  num_gaussians = numel(locations)/2;

  % Extract the locations
  locations = round(crsUnitToUnit(locations_units,locations,CRS.PIXELUNIT));
  if     all( size(locations) == [num_gaussians,2] )
    mean_1 = locations(:,1);
    mean_2 = locations(:,2);
  elseif all( size(locations) == [2,num_gaussians] )
    mean_1 = locations(1,:);
    mean_2 = locations(2,:);
  else
    error('locations must be an array of length 2 vectors, and there must be the same number of locations as there are values in the other vectorised parameters.');
  end

  % Extract the deviations
  deviations = round(crsUnitToUnit(deviations_units,deviations,CRS.PIXELUNIT));
  if      all( size(deviations) == [num_gaussians,2] )
    standard_deviation_1 = deviations(:,1);
    standard_deviation_2 = deviations(:,2);
  elseif  all( size(deviations) == [2,num_gaussians] )
    standard_deviation_1 = deviations(1,:);
    standard_deviation_2 = deviations(2,:);
  else
    error('deviations must be an array of length 2 vectors, and there must be the same number of deviations as there are values in the other vectorised parameters.');
  end

  % Extract the correlations
  if min(size(correlations)) ~= 1 || max(size(correlations)) ~= num_gaussians
    error('correlations must be a vector, and there must be the same number of correllations as there are values in the other vectorised parameters.');
  end
  cross_correlation_coefficient = correlations;
  if min(cross_correlation_coefficient) < 0 || 1 <= max(cross_correlation_coefficient)
    error('The cross correlation coefficients should lie between 0 and 1: 0 <= rho < 1.');
  end

  % Extract the angles
  if min(size(angles)) ~= 1 || max(size(angles)) ~= num_gaussians
    error('angles must be a vector, and there must be the same number of angles as there are values in the other vectorised parameters.');
  end

  % Extract the scaling_factors
  if min(size(scaling_factors)) ~= 1 || max(size(scaling_factors)) ~= num_gaussians
    error('scaling_factors must be a vector, and there must be the same number of scaling_factors as there are values in the other vectorised parameters.');
  end

  % Unfortunately, we need a loop to add all our gaussians together, as
  % vectorising it would require that we held all of our gaussians in
  % memory at once - not really a practical proposition.
  gaussian_mesh = zeros(size(x1_grid));
  for i = 1:num_gaussians

    % Break the angle down into two orthogonal components;
    [rotation_factor_11,rotation_factor_12] = pol2cart( angles(i),      1);
    [rotation_factor_21,rotation_factor_22] = pol2cart((angles(i))+pi/2,1);

    % Rotate the coordinate grid about the center of the gaussian
    % ( this also removes the need to talk about the mean in the
    %   function - it will always be zero! )
    x1 =   ( rotation_factor_11 .* (x1_grid - mean_1(i)) ) ...
         + ( rotation_factor_12 .* (x2_grid - mean_2(i)) ) ;

    x2 =   ( rotation_factor_21 .* (x1_grid - mean_1(i)) ) ...
         + ( rotation_factor_22 .* (x2_grid - mean_2(i)) ) ;

    % Call our parameters something simple for brevety's sake.
    sd1 = standard_deviation_1(i);
    sd2 = standard_deviation_2(i);

    rho = cross_correlation_coefficient(i);
    k   = scaling_factors(i);


    % -------------------------------------------------------------------------
    % We could use the following formula to calculate a 'proper gaussian'
    % (where the area under the surface sums to 1):
    %
    % z1  = ((x1 - m1).^2) / sd1 .^ 2;
    % z2  = ((x2 - m2).^2) / sd2 .^ 2;
    % z3  = ( 2 .* rho .* (x1 - m1) .* (x2 - m2) ) / (sd1 .* sd2);
    % z   = z1 + z2 - z3;
    % exponent = - (   z / (  2 .* ( 1 - (rho.^2) )  )   );
    % constant = 1 / (  2 * pi * sd1 * sd2 * sqrt( 1 - (rho.^2) ) );
    % gaussian_mesh = constant .* exp(exponent);
    %
    % gaussian_mesh = gaussian_mesh + (1/(2*pi*sd1*sd2*sqrt(1-(rho.^2)))).*exp(-(((((x1-m1).^2)/sd1.^2)+(((x2-m2).^2)/sd2.^2)-((2.*rho.*(x1-m1).*(x2-m2))/(sd1.*sd2)))/(2.*(1-(rho.^2)))));

    % But we are only generally interested in the shape, so we should
    % probably use instead the simpler formula, where each gaussian
    % sums to 1:
    %
    % z1  = ((x1 - m1).^2) / sd1 .^ 2;
    % z2  = ((x2 - m2).^2) / sd2 .^ 2;
    % z3  = ( 2 .* rho .* (x1 - m1) .* (x2 - m2) ) / (sd1 .* sd2);
    % z   = z1 + z2 - z3;
    % gaussian_mesh = exp(-z);
    %
    % gaussian_mesh = gaussian_mesh + exp(( (((x1 - m1).^2) / sd1 .^ 2) + (((x2 - m2).^2) / sd2 .^ 2) - (( 2 .* rho .* (x1 - m1) .* (x2 - m2) ) / (sd1 .* sd2)) ))
    %
    % We can simplify this further by dealing with the means when we deal
    % with coordinate system rotation, leaving us with:
    % -------------------------------------------------------------------------

    gaussian_mesh = gaussian_mesh + ( k * exp(                           ...
                                               -(   ((x1.^2)/sd1.^2)     ...
                                                  + ((x2.^2)/sd2.^2)     ...
                                                  - (   (2.*rho.*x1.*x2) ...
                                                      / (sd1.*sd2)       ))));

  end
