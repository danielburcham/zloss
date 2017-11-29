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