function [zloss,Ad] = zloss(geoimage,tomogram,varargin)
%{
DESCRIPTION
Using PiCUS sonic tomograms, this function computes the percent loss to
section modulus, Z_{LOSS} (%), for trees affected by internal defects and 
visualizes the associated estimates.

NOTES
This version computes estimates for trees with one or more damaged parts
and one or more open cavities. When the spatial extent of pixels is not 
provided as an input parameter, the two input images must be the same size. 

[figuretype,zloss,Ad] = z_loss(geoimage,tomogram);
[figuretype,zloss,Ad] = z_loss(geoimage,tomogram,'parameter1',value1,...);

REQUIRED INPUTS
geoimage: string - filepath for geometry image file showing blue trunk
boundary line
tomogram: string - filepath for tomogram showing visualized damage pattern

OPTIONAL INPUTS
colors: string - colors used to select damaged parts, either blue, violet,
and green ('BVG') or blue and violet ('BV').
geopixelextentinworldx: float - spatial extent of each pixel in the
x-direction for the geometry image file
geopixelextentinworldy: float - spatial extent of each pixel in the
y-direction for the geometry image file
tomopixelextentinworldx: float - spatial extent of each pixel in the
x-direction for the tomogram
tomopixelextentinworldy: float - spatial extent of each pixel in the
y-direction for the tomogram

PARAMETER                     CLASS       DEFAULT VALUE
---------------------------------------------------------------------------
colors                        string      'BVG'
geopixelextentinworldx        float       1
geopixelextentinworldy        float       1
tomopixelextentinworldx       float       1
tomopixelextentinworldy       float       1

OUTPUTS
zloss: mx2 matrix containing estimates of degrees clockwise rotation from
vertical and the percent loss to Z, ZLOSS (%) in the first and second 
column, respectively.
Ad: Percent (%) of total damaged cross sectional area.

---------------------------------------------------------------------------
%}
format long

% Check function call
if nargin < 2
    error('z_loss requires at least two inputs');
elseif ~isempty(varargin)
    if mod(length(varargin),2)~=0
        if (length(varargin)==1 && isstruct(varargin{1}))
            varargin = reshape([fieldnames(varargin{1})...
                struct2cell(varargin{1})]',1,[]);
        else
            error(strcat('Inputs must be paired: z_loss(geoimage,',...
                'tomogram,''PropertyName'',PropertyValue,...)'));
        end
    elseif ~any(strcmp(varargin{find(strcmp(varargin,'colors'))+1},...
            {'BVG','BV'}))
        error('colors must be set equal to ''BV'' or ''BVG''');
    end
end

% Default parameters
Colors = 'BVG';
pixelExtentInWorld11 = 1;
pixelExtentInWorld12 = 1;
pixelExtentInWorld21 = 1;
pixelExtentInWorld22 = 1;

% User-specified parameters
if ~isempty(varargin)
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'colors'
                Colors = varargin{i+1};
            case 'geopixelextentinworldx'
                pixelExtentInWorld11 = varargin{i+1};
            case 'geopixelextentinworldy'
                pixelExtentInWorld12 = varargin{i+1};
            case 'tomopixelextentinworldx'
                pixelExtentInWorld21 = varargin{i+1};
            case 'tomopixelextentinworldy'
                pixelExtentInWorld22 = varargin{i+1};
        otherwise
            error([varargin{i} 'is not a valid property for the z_loss function.']);
        end
    end
end
input1 = find(strcmp({'BV','BVG'},Colors))-1;

% Import tomogram image file (.jpg)
A = imread(geoimage,'jpg');
B = imread(tomogram,'jpg');
if all([pixelExtentInWorld11,pixelExtentInWorld12,pixelExtentInWorld21,...
        pixelExtentInWorld22]==1) && ~isequal(size(A),size(B))
    error(strcat('Either the two images must be the same size or ',...
        'the pixel extents must be specified.'));
end
sz1 = size(A);
R1 = imref2d(sz1,pixelExtentInWorld11,pixelExtentInWorld12);
sz2 = size(B);
R2 = imref2d(sz2,pixelExtentInWorld21,pixelExtentInWorld22);

% Segment image
O = segment(A,B,R1,R2,input1);

% Search for open cavities and modify outer boundary
for i=1:size(O{2,1},1)-1
    [in,~]=inpolygon(O{2,1}{i+1,1}(:,1),O{2,1}{i+1,1}(:,2),...
        O{1,1}(:,1),O{1,1}(:,2));
    if any(~in)
        [x,y] = polybool('-',O{2,1}{1,1}(:,1),O{2,1}{1,1}(:,2),...
            O{2,1}{i+1,1}(:,1),O{2,1}{i+1,1}(:,2));
        O{2,1}{1,1} = [x,y];
        O{2,1}{i+1,1}=[];
        clear x y
    end
end
if isShapeMultipart(O{2,1}{1,1}(:,1),O{2,1}{1,1}(:,2))
    [X,Y]=polysplit(O{2,1}{1,1}(:,1),O{2,1}{1,1}(:,2));
    U=cell(sum(cellfun('length',X)>=75),1);
    ia=1;
    for i=1:length(X)
        if length(X{i,1}) > 75
            U{ia,1}=[X{i,1} Y{i,1}];
            ia=ia+1;
        end
    end
    O{2,1}{1,1}=U;
else
    U=cell(1,1);
    U{1,1} = O{2,1}{1,1};
    O{2,1}{1,1}=U;
end

% Reshape O to remove empty cells and do isShape
O{2,1} = O{2,1}(~cellfun('isempty',O{2,1}));

% Compute section properties and percent reduction to Z
[zloss,Ad,L] = mcl(O);

% Visualize results
if ~iscell(O{2,1}{1,1})
    M = zeros(length(zloss(:,2)),3);
    maxDistance = zeros(1,length(O{1,1}(:,1)));
    for k = 1:length(O{1,1}(:,1))
        distances = sqrt((O{1,1}(:,1)-O{1,1}(k,1)).^2+(O{1,1}(:,2)-O{1,1}(k,2)).^2);
        maxDistance(k) = max(distances);
    end
    rad = 0.65*max(maxDistance);
    ang = zloss(3,1)-zloss(2,1);
    theta = pi/2:-ang:(-3*pi/2)+ang;
    M(:,1) = rad.*cos(theta)+L(1);
    M(:,2) = rad.*sin(theta)+L(2);
    M(:,3) = zloss(:,2);
    patch(M(:,1),M(:,2),M(:,3),'facecolor','none','edgecolor','interp',...
        'linewidth',2);
    colormap('jet');
    hold on;
    for ia=1:size(O{2,1},1)
        plot(O{2,1}{ia,1}(:,1),O{2,1}{ia,1}(:,2),'Color',[0 0 0],'LineWidth',0.75);
    end
    c = colorbar;
    c.Label.String = 'Percent reduction to section modulus, \itZ_{LOSS}\rm (%)';
    c.Label.FontSize = 12;
elseif iscell(O{2,1}{1,1})
    M = zeros(length(zloss(:,2)),3);
    maxDistance = zeros(1,length(O{1,1}(:,1)));
    for k = 1:length(O{1,1}(:,1))
        distances = sqrt((O{1,1}(:,1)-O{1,1}(k,1)).^2+(O{1,1}(:,2)-O{1,1}(k,2)).^2);
        maxDistance(k) = max(distances);
    end
    rad = 0.65*max(maxDistance);
    ang = zloss(3,1)-zloss(2,1);
    theta = pi/2:-ang:(-3*pi/2)+ang;
    M(:,1) = rad.*cos(theta)+L(1);
    M(:,2) = rad.*sin(theta)+L(2);
    M(:,3) = zloss(:,2);
    patch(M(:,1),M(:,2),M(:,3),'facecolor','none','edgecolor','interp',...
        'linewidth',2);
    colormap('jet');
    hold on;
    for ia=1:size(O{2,1},1)-1+size(O{2,1}{1,1},1)
        if ia <= size(O{2,1}{1,1},1)
            plot(O{2,1}{1,1}{ia,1}(:,1),O{2,1}{1,1}{ia,1}(:,2),...
                'Color',[0 0 0],'LineWidth',0.75);
        elseif ia > size(O{2,1}{1,1},1)
            plot(O{2,1}{ia-size(O{2,1}{1,1},1)+1,1}(:,1),...
                O{2,1}{ia-size(O{2,1}{1,1},1)+1,1}(:,2),'Color',...
                [0 0 0],'LineWidth',0.75);
        end
    end
    c = colorbar;
    c.Label.String = 'Percent reduction to section modulus, \itZ_{LOSS}\rm (%)';
    c.Label.FontSize = 12;
end

end

function [N] = segment(A,B,R1,R2,input1)
%{
DESCRIPTION
Segments tomogram into solid and decayed regions and returns the associated 
boundary coordinates

NOTES
This version computes estimates for trees with one more more columns of 
decay and one or more open cavities

---------------------------------------------------------------------------
INPUTS
A: Geometry image file showing only the blue trunk boundary line.
B: Sonic tomogram image file showing visualized decay pattern.
R1: Coordinate systen reference object for trunk geometry file.
R2: Coordinate systen reference object for tomogram file.
input1: Choose liberal (0) or conservative (1) estimate.

OUTPUTS
N: 1x2 cell array with Cartesian coordinates (m) for solid and decayed
regions contained in the first and second cell, respectively. Coordinate
pairs are recorded as row-wise entries in a mX2 matrix.
---------------------------------------------------------------------------
%}

% Choose histogram thresholds for solid region
C = zeros(6,1);
C(1,1) = 0.014;
C(2,1) = 0.776;
C(3,1) = 0.160;
C(4,1) = 1.000;
C(5,1) = 0.000;
C(6,1) = 1.000;

% Convert RGB image to HSV color space
D = rgb2hsv(A);

% Create mask based on histogram thresholds for solid region
E = (D(:,:,1) >= C(1,1) ) & (D(:,:,1) <= C(2,1)) & ...
    (D(:,:,2) >= C(3,1) ) & (D(:,:,2) <= C(4,1)) & ...
    (D(:,:,3) >= C(5,1) ) & (D(:,:,3) <= C(6,1));

% Find connected regions, select maximum, and trace boundary
F = bwconncomp(E,4);
numPixels = cellfun(@numel,F.PixelIdxList);
[~,idx] = max(numPixels);
G = zeros(size(A,1),size(A,2));
G(F.PixelIdxList{idx}) = 1;
[H,~,~,I] = bwboundaries(G,'holes',8);

%Convert RGB image to LAB color space
Z = rgb2lab(B);

% Choose histogram thresholds for decayed (blue and purple) region
Y = zeros(6,1);
Y(1,1) = 20.922;
Y(2,1) = 100.000;
Y(3,1) = -26.902;
Y(4,1) = 81.975;
Y(5,1) = -63.537;
Y(6,1) = -3.780;
X = (Z(:,:,1) >= Y(1,1) ) & (Z(:,:,1) <= Y(2,1)) & ...
    (Z(:,:,2) >= Y(3,1) ) & (Z(:,:,2) <= Y(4,1)) & ...
    (Z(:,:,3) >= Y(5,1) ) & (Z(:,:,3) <= Y(6,1));
W = bwareaopen(X,30);

if input1==1
    % Choose histogram thresholds for decayed (green) region
    V=zeros(6,1);
    V(1,1)=0.000;
    V(2,1)=100.000;
    V(3,1)=-68.602;
    V(4,1)=-10.461;
    V(5,1)=-99.915;
    V(6,1)=62.944;
    U = (Z(:,:,1) >= V(1,1) ) & (Z(:,:,1) <= V(2,1)) & ...
        (Z(:,:,2) >= V(3,1) ) & (Z(:,:,2) <= V(4,1)) & ...
        (Z(:,:,3) >= V(5,1) ) & (Z(:,:,3) <= V(6,1));
    T = bwareaopen(U,30);
    S = logical(W+T);
    W = S;
end

R = bwconncomp(W,6);
numPixels = cellfun(@numel,R.PixelIdxList);
Q = find(numPixels>100); % Select decayed areas with > 100 pixels
P = zeros(size(B,1),size(B,2));
for i = 1:size(Q,2)
    P(R.PixelIdxList{Q(i)}) = 1;
end
P(1:20,:) = 0; % Exclude tomogram colorbar from decayed areas
P=logical(P);
O = bwboundaries(P,'noholes',8);

% Convert to registered Cartesian coordinate system
ind1 = find(I(:,1));
[~,ind2] = max(cellfun('length',H(find(I(:,1)))));
ind3 = ind1(ind2);
[x1, y1] = intrinsicToWorld(R1,H{ind3,1}(:,2),H{ind3,1}(:,1)); %Solid region
ymax = max(y1);
y1 = ymax - y1;
N = cell(2,1);
N{1,1} = [x1, y1];
N{2,1}{1,1} = [x1, y1];
if size(O,1)==1
    [x2, y2] = intrinsicToWorld(R2,O{1,1}(:,2),O{1,1}(:,1));
    y2 = ymax - y2;
    N{2,1}{2,1} = [x2, y2];
elseif size(O,1)>1
    for i=1:size(O,1)
        [x2, y2] = intrinsicToWorld(R2,O{i,1}(:,2),O{i,1}(:,1));
        y2 = ymax - y2;
        N{2,1}{i+1,1} = [x2, y2];
    end
end
end

function [H,I,L] = mcl(X)
%{
Description:
Computes percent loss to section modulus, Z (m^3), for trees affected by 
decay

Notes:
This version computes estimates for trees with one more more columns of 
decay and one or more open cavities

---------------------------------------------------------------------------
Inputs:
X: 1x2 cell array with Cartesian coordinates (m) for solid and decayed 
sections contained in the first and second cell, respectively. Coordinate 
pairs are recorded as row-wise entries in a mX2 matrix.

Outputs:
H: mx2 matrix containing estimates of the percent loss to Z and degrees 
clockwise rotation from North in the first and second column, respectively.
I: Percent (%) of cross-sectional area occupied by decay.
L: Centroid coordinates for hollow section.

---------------------------------------------------------------------------
%}

% Pre-allocate matrices
ns = size(X{2,1}{1,1},1); % Number of solid shapes comprising hollow section
nv = size(X{2,1},1)-1; % Number of void shapes comprising hollow section
z_loss=zeros(length(X{1,1}(:,1)),1);
ang = 2*pi/length(X{1,1}(:,1));
B=zeros(length(X{1,1}(:,1)),ns+nv+1); % Inertial moments: all component shapes
C=zeros(length(X{1,1}(:,1)),ns+nv+1); % Areas: all component shapes
D=zeros(length(X{1,1}(:,1)),2*(ns+nv+1)); % Centroid coorindates: all component shapes
E=zeros(length(X{1,1}(:,1)),2); % Centroid coordinates: hollow section
F=zeros(length(X{1,1}(:,1)),2*(ns+nv)); % Area moments in original coordinates
G=zeros(length(X{1,1}(:,1)),1); % Inertial moments: hollow section
J=zeros(length(X{1,1}(:,1)),1); % y: hollow section
K=zeros(length(X{1,1}(:,1)),1); % y: solid section

for i=1:length(X{1,1}(:,1))
    % Determine I, A, and centroid of all component shapes
    for ia = 1:ns+nv+1
        if ia == 1
            x = X{1,1}(:,1);
            y = X{1,1}(:,2);
        elseif ia >= 2 && ia <= ns+1
            x = X{2,1}{1,1}{ia-1,1}(:,1);
            y = X{2,1}{1,1}{ia-1,1}(:,2);
        elseif ia > ns+1
            x = X{2,1}{ia-ns,1}(:,1);
            y = X{2,1}{ia-ns,1}(:,2);
        end
        if ~ispolycw(x,y)
            [x,y]=poly2cw(x,y);
        end
        
        % Number of vertices
        [x,~] = shiftdim(x);
        [y,~] = shiftdim(y);
        [n,~] = size(x);
        
        % Temporarily shift data to mean of vertices for improved accuracy
        xm = mean(x);
        ym = mean(y);
        x = x - xm*ones(n,1);
        y = y - ym*ones(n,1);
        
        % Delta x and delta y
        dx = x( [ 2:n 1 ] ) - x;
        dy = y( [ 2:n 1 ] ) - y;
        
        % Summations for CW boundary integrals
        A = sum( y.*dx - x.*dy )/2;
        Axc = sum( 6*x.*y.*dx -3*x.*x.*dy +3*y.*dx.*dx +dx.*dx.*dy )/12;
        Ayc = sum( 3*y.*y.*dx -6*x.*y.*dy -3*x.*dy.*dy -dx.*dy.*dy )/12;
        Ixx = sum( 2*y.*y.*y.*dx -6*x.*y.*y.*dy -6*x.*y.*dy.*dy ...
            -2*x.*dy.*dy.*dy -2*y.*dx.*dy.*dy -dx.*dy.*dy.*dy )/12;
        
        % Centroidal moments
        xc = Axc / A;
        yc = Ayc / A;
        Iuu = Ixx - A*yc*yc;
        
        % Replace mean of vertices
        x_cen = xc + xm;
        y_cen = yc + ym;
        
        B(i,ia) = Iuu;
        C(i,ia) = A;
        D(i,2*ia-1) = x_cen;
        D(i,2*ia) = y_cen;
    end
    
    B(i,ns+2:end) = -1.*B(i,ns+2:end);
    C(i,ns+2:end) = -1.*C(i,ns+2:end);
    
    % Determine centroid of hollow section
    O = vertcat(X{2,1}{1,1},X{2,1}(2:end,1));
    M = cellfun(@sum,O,'UniformOutput',false);
    P = cellfun(@length,O,'UniformOutput',false);
    xm = sum(cellfun(@(c) c(1,1), M(:,1)))/sum(cellfun(@sum,P));
    ym = sum(cellfun(@(c) c(1,2), M(:,1)))/sum(cellfun(@sum,P));
    for ia = 1:ns+nv
        if ia <= ns
            x = X{2,1}{1,1}{ia,1}(:,1);
            y = X{2,1}{1,1}{ia,1}(:,2);
        elseif ia > ns
            x = X{2,1}{ia-ns+1,1}(:,1);
            y = X{2,1}{ia-ns+1,1}(:,2);
        end
        
        if ~ispolycw(x,y)
            [x,y]=poly2cw(x,y);
        end
        
        % Number of vertices
        [x,~] = shiftdim(x);
        [y,~] = shiftdim(y);
        [n,~] = size(x);
        
        % Temporarily shift data to mean of vertices for improved accuracy
        x = x - xm*ones(n,1);
        y = y - ym*ones(n,1);
        
        % Delta x and delta y
        dx = x( [ 2:n 1 ] ) - x;
        dy = y( [ 2:n 1 ] ) - y;
        
        % Summations for CW boundary integrals
        Axc = sum( 6*x.*y.*dx -3*x.*x.*dy +3*y.*dx.*dx +dx.*dx.*dy )/12;
        Ayc = sum( 3*y.*y.*dx -6*x.*y.*dy -3*x.*dy.*dy -dx.*dy.*dy )/12;
        
        F(i,2*ia-1) = Axc;
        F(i,2*ia) = Ayc;
    end
    
    F(i,2*ns+1:end) = -1.*F(i,2*ns+1:end);
    
    E(i,1)=sum(F(i,1:2:end))./sum(C(i,2:end)) + xm; % Xbar
    E(i,2)=sum(F(i,2:2:end))./sum(C(i,2:end)) + ym; % Ybar
    
    %Calculate Iuu for hollow section using parallel axis theorem
    G(i) = sum(B(i,2:end)+(C(i,2:end).*(D(i,4:2:end)-E(i,2)).^2));
    % where the first and second squared terms are the vertical distance 
    % between the axis of rotation and the centroid of the solid and
    % decayed area, respectively
    
    % Determine maximum vertical distance between hollow section neutral 
    % axis and upper boundary
    u = size(X{2,1}{1,1},1);
    for ia=1:size(X{2,1}{1,1},1)
        u(ia)=max(abs(X{2,1}{1,1}{ia,1}(:,2)-E(i,2)));
    end
    J(i) = max(u);
    
    % Determine maximum vertical distance between solid section neutral 
    % axis and upper boundary
    K(i)=max(abs(X{1,1}(:,2)-D(i,2)));
    
    % Calculate loss of section modulus and store in array
    z_loss(i) = ((B(i,1)./K(i))-(G(i)./J(i)))./(B(i,1)./K(i));
    
    % Rotate coordinate system by 1 step about the centroid of the hollow
    % section
    for ie=1:ns+nv+1
        if ie == 1
            r = [X{1,1}(:,1) X{1,1}(:,2)]';
            R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
            v = repmat([D(1,1); D(1,2)], 1, size(r,2));
            y = (R*(r-v)+v)';
            X{1,1}(:,1) = y(:,1);
            X{1,1}(:,2) = y(:,2);
        elseif ie >= 2 && ie <= ns+1
            r = [X{2,1}{1,1}{ie-1,1}(:,1) X{2,1}{1,1}{ie-1,1}(:,2)]';
            R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
            v = repmat([E(1,1); E(1,2)], 1, size(r,2));
            y = (R*(r-v)+v)';
            X{2,1}{1,1}{ie-1,1}(:,1) = y(:,1);
            X{2,1}{1,1}{ie-1,1}(:,2) = y(:,2);
        elseif ie > ns+1
            r = [X{2,1}{ie-ns,1}(:,1) X{2,1}{ie-ns,1}(:,2)]';
            R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
            v = repmat([E(1,1); E(1,2)], 1, size(r,2));
            y = (R*(r-v)+v)';
            X{2,1}{ie-ns,1}(:,1) = y(:,1);
            X{2,1}{ie-ns,1}(:,2) = y(:,2);
        end
    end
    
end

% Output
theta = 0:ang:(2*pi)-ang;
theta = theta';
H = [theta, z_loss];
I = 1-sum(mean(C(:,2:end)))/mean(C(:,1));
L = [D(1,1), D(1,2)];

end