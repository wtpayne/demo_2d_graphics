function IMG = crsCreateStringPrimitive(STRING,SIZE,FONT)
% CRSCREATESTRINGPRIMITIVE creates an image matrix of a specified string
% -----------------------------------------------------------------------------
% crsCreateStringPrimitive
% ========================
% This function creates an image matrix (values in the range 0,1) from a 
% specified string. The PointSize parameter specifies the size in points
% (i.e. 12 for a 12-point font).
%
% Note: This function works by capturing the string from a MATLAB window
%       using 'getframe'. It therefore has limitations on the length of 
%       string that it can accept (related to your windows desktop resolution)
%       and also momentarily creates a MATLAB figure window to capture the
%       string.
%
%       As an alternative, the crsCreateGlyphPrimitive function gets a single 
%       character image matrix directly from the Windows API.
%
% Usage
% =====
% ImageMatrix = crsCreateStringPrimitive(String);
% ImageMatrix = crsCreateStringPrimitive(String,PointSize);
% ImageMatrix = crsCreateStringPrimitive(String,PointSize,FontName);
%
% Reference page in Help browser
% <a href="matlab:web(['jar:file:',which('crsHelp\help.jar'),'!/crs\tools\stimulusmatrices\crsCreateStringPrimitive.html'],'-helpbrowser');">crsCreateStringPrimitive HTML help.</a>                                                                                                                       
%
% -----------------------------------------------------------------------------

      if nargin == 1
    SIZE = 40;
    FONT = 'Helvetica';
  elseif nargin==2  
    FONT = 'Helvetica';
  end  

% Get some information about the root object.
TextBoxWidth  = 400;
TextBoxHeight = 80; 
ScreenSize = get(0,'ScreenSize');
TextBoxX = ScreenSize(3) - TextBoxWidth;
TextBoxY = ScreenSize(4) - TextBoxHeight;

% Make a figure containing some axes.
FIG = figure('Visible','off','Name','String Capture Window.','Units','pixels', ...
             'ToolBar','none','MenuBar','none','Position',[TextBoxX,TextBoxY,TextBoxWidth,TextBoxHeight]);
AXS = axes('Units','pixels','Position',[0,0,TextBoxWidth,TextBoxHeight],'Color',[1,1,1]);

% Create some text
TXT    = text(50,50,STRING,'Units','pixels',            ...
                           'FontAngle'  , 'normal',     ...
                           'FontName'   , FONT,         ...
                           'FontSize'   , SIZE,         ...
                           'FontUnits'  , 'points',     ...
                           'FontWeight' , 'normal'); 
TextExtent = get(TXT,'Extent');
TextBoxWidth  = TextExtent(3) + TextExtent(1);
TextBoxHeight = TextExtent(2) + TextExtent(4);
ScreenSize = get(0,'ScreenSize');
TextBoxX = ScreenSize(3) - TextBoxWidth;
TextBoxY = ScreenSize(4) - TextBoxHeight;

% Make a figure containing some axes.
set(FIG,'Position',[TextBoxX,TextBoxY,TextBoxWidth,TextBoxHeight]);
set(AXS,'Position',[0,0,TextBoxWidth,TextBoxHeight]);

% Get the image from the frame.
         padding = [-10,-10,20,20];
IMG    = getframe(FIG,get(TXT,'Extent')+padding); 
delete(FIG);drawnow;
IMG    = double(IMG.cdata(:,:,1));

% Scale image to between 0 and 1.
IMG    = IMG-min(IMG(:));
IMG    = IMG/max(IMG(:));

% Invert so background is 0 and text is 1.
IMG    = 1-IMG;
