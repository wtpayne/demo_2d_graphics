function varargout = utilCRSgenCoordsCart(varargin)
% This function creates a cartesian coordinate system of a specified size.
% -----------------------------------------------------------------------------
% utilCRSgenCoordsCart
% ====================
%  Often, when creating visual stimuli, we want to express our stimulus as a
%  function of spatial position. i.e; in the following terms:
%
%    stim = f(x,y).
%
%  In order to do this in MATLAB, we need to create two matrices [x,y] that
%  are the same size as our desired stimulus matrix. Each point in our stimulus
%  matrix will therefore have two corresponding points, x and y.
%
%  Usage:
%  ======
%   If no input parameters are provided, this function assumes that you want
%   a coordinate-system the size of the stimulus display. (units in pixels)
%
%    [x,y] = crsGenCoordsCart;
%    [x,y] = crsGenCoordsCart([width,height]);
%    [x,y] = crsGenCoordsCart( width,height );
%
% Reference page in Help browser
% <a href="matlab:web(['jar:file:',which('crsHelp\help.jar'),'!/crs\tools\stimulusmatrices\coordinates\utilCRSgenCoordsCart.html'],'-helpbrowser');">utilCRSgenCoordsCart HTML help.</a>                                                                                                                   
%
% -----------------------------------------------------------------------------

  % First of all, we need to extract the width and the height
  % from the input parameters.
  if nargin == 0
    width  = crsGetScreenWidthPixels;
    height = crsGetScreenHeightPixels;
  elseif nargin == 1 && numel(varargin{1}) == 2
    sizev = varargin{1};
    width  = sizev(1);
    height = sizev(2);
  elseif nargin == 2 && crsIsScalar(varargin{1}) && crsIsScalar(varargin{2})
    width  = varargin{1};
    height = varargin{2};
  else
    error('Unrecognised input parameters');
  end

  % Now, we create a cartesian coordinate system with the origin at the
  % centre of the screen. Each point on the screen needs two values: an
  % X-coordinate and a Y-coordinate; This is provided by two matrices of
  % equal size.
  halfX = (width / 2);
  halfY = (height / 2);
  xVector = linspace(-halfX,halfX,width);
  yVector = linspace(-halfY,halfY,height);
  % We need to flip y upside down because matlab is weird.
  yVector = fliplr(yVector);
  [xGrid,yGrid] = meshgrid(xVector,yVector);

  % Now, we format the output parameters as they have been requested.
  if nargout == 0
  elseif nargout == 1
    coords.X = xGrid;
    coords.Y = yGrid;
    varargout{1} = coords;
  elseif nargout == 2
    varargout{1} = xGrid;
    varargout{2} = yGrid;
  else
    error('Unrecognised output parameters');
  end
