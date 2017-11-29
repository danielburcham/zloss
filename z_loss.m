%{
Description:
Using PiCUS sonic tomograms, this script computes the percent loss to
section modulus, Z (m^3), for trees affected by decay and visualizes the
associated estimates.

Notes:
This version computes estimates for trees with one or more columns of decay
and one or more open cavities. For sonic tomograms, liberal or conservative
estimates correspond to selecting decayed areas by excluding or including 
green, respectively.

---------------------------------------------------------------------------
Inputs:
input1: Choose liberal (0) or conservative (1) estimate.
A: Geometry image file showing only the blue trunk boundary line.
B: Sonic tomogram image file showing visualized decay pattern.
pixelExtentInWorld11: Spatial extent of each pixel in the x-direction for
the geometry image file.
pixelExtentInWorld12: Spatial extent of each pixel in the y-direction for
the geometry image file.
pixelExtentInWorld21: Spatial extent of each pixel in the x-direction for
the tomogram image file.
pixelExtentInWorld22: Spatial extent of each pixel in the y-direction for
the tomogram image file.

Outputs:
H: mx2 matrix containing estimates of degrees clockwise rotation from
vertical and the percent loss to Z in the first and second column,
respectively.
I: Percent (%) of cross sectional area occupied by decay.

---------------------------------------------------------------------------
%}
format long

% Import tomogram image file (.jpg)
input1 = input('Choose liberal (0) or conservative (1) estimate: ');
filename1 = input('Specify geometry file path: ','s');
A = imread(filename1,'jpg');
sz1 = size(A);
pixelExtentInWorld11 = input('Specify pixel extent in x-direction (m): ');
pixelExtentInWorld12 = input('Specify pixel extent in y-direction (m): ');
R1 = imref2d(sz1,pixelExtentInWorld11,pixelExtentInWorld12);
filename2=input('Specify tomogram file path: ','s');
B = imread(filename2,'jpg');
sz2 = size(B);
pixelExtentInWorld21 = input('Specify pixel extent in x-direction (m): ');
pixelExtentInWorld22 = input('Specify pixel extent in y-direction (m): ');
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
[H,I,L] = mcl(O);

% Visualize results
if ~iscell(O{2,1}{1,1})
    M = zeros(length(H(:,2)),3);
    maxDistance = zeros(1,length(O{1,1}(:,1)));
    for k = 1:length(O{1,1}(:,1))
        distances = sqrt((O{1,1}(:,1)-O{1,1}(k,1)).^2+(O{1,1}(:,2)-O{1,1}(k,2)).^2);
        maxDistance(k) = max(distances);
    end
    rad = 0.65*max(maxDistance);
    ang = H(3,1)-H(2,1);
    theta = pi/2:-ang:(-3*pi/2)+ang;
    M(:,1) = rad.*cos(theta)+L(1);
    M(:,2) = rad.*sin(theta)+L(2);
    M(:,3) = H(:,2);
    patch(M(:,1),M(:,2),M(:,3),'facecolor','none','edgecolor','interp',...
        'linewidth',2);
    colormap('jet');
    hold on;
    for ia=1:size(O{2,1},1)
        plot(O{2,1}{ia,1}(:,1),O{2,1}{ia,1}(:,2),'Color',[0 0 0],'LineWidth',0.75);
    end
    c = colorbar;
    c.Label.String = 'Percent reduction to section modulus, \itZ\rm (m^{3})';
    c.Label.FontSize = 12;
elseif iscell(O{2,1}{1,1})
    M = zeros(length(H(:,2)),3);
    maxDistance = zeros(1,length(O{1,1}(:,1)));
    for k = 1:length(O{1,1}(:,1))
        distances = sqrt((O{1,1}(:,1)-O{1,1}(k,1)).^2+(O{1,1}(:,2)-O{1,1}(k,2)).^2);
        maxDistance(k) = max(distances);
    end
    rad = 0.65*max(maxDistance);
    ang = H(3,1)-H(2,1);
    theta = pi/2:-ang:(-3*pi/2)+ang;
    M(:,1) = rad.*cos(theta)+L(1);
    M(:,2) = rad.*sin(theta)+L(2);
    M(:,3) = H(:,2);
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
    c.Label.String = 'Percent reduction to section modulus, \itZ\rm (m^{3})';
    c.Label.FontSize = 12;
end
