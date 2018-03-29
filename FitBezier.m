function pointsOut = FitBezier(pointsIn)
    maxError = 20;%500;
    nloops = size(pointsIn,1);
    pointsOut = [];
    
    % for each loop
    for ii = 1:nloops
        % select iith loop
        points = pointsIn{ii};
        % split curve in half
        nrows = size(points,1);
        points1 = points(1:floor(nrows/2 + 1),:);
        points2 = points(floor(nrows/2 + 1):end,:);
        
        % for each half
        for jj = 1:2
            if jj == 1
                bpoints = points1;
            else
                bpoints = points2;
            end
            % create empty array for Bezier points
            pBez = zeros(4,2);
            % take endpoints as p0 and p3
            pBez(1,:) = bpoints(1,:);
            pBez(2,:) = bpoints(1,:);
            pBez(3,:) = bpoints(end,:);
            pBez(4,:) = bpoints(end,:);
            % function to be minimized
            error = @(c) sumSquaredError(findErrors(bpoints, [pBez(1,:);c(1,:);c(2,:);pBez(4,:)]));
            % minimize error
            c0 = [pBez(2,:);pBez(3,:)];
            options = optimoptions('fminunc','Algorithm','quasi-newton');
            [cpoints, e] = fminunc(error,c0,options);
            if e < maxError
                pBez(2:3,:) = cpoints;
                pointsOut = cat(3,pointsOut,pBez);
            else
                newInput = cell(1);
                newInput{1} = bpoints;
                pointsOut = cat(3,pointsOut,FitBezier(newInput));
            end
        end % jj
    end % ii
end % main

% Evaluates the locaion of a point and the tangent vector on a cubic Bezier
% curve
% INPUT: 
%   pBez - a 4x3 matrix containing the points b0, b1, b2, and b3 as row
%   vectors
%   u - parameter value
%OUTPUT:
%   p - point location on curve
%   v - tangent vector to the curve
function [p, v] = bezier(pBez, u)
    u = [u^3 u^2 u 1];

    M = [-1  3 -3  1;
          3 -6  3  0;
         -3  3  0  0;
          1  0  0  0];
      
    D = [ 0  0  0  0;
         -3  9 -9  3;
          6 -12 6  0;
         -3  3  0  0];
     
    p = u*M*pBez;
    v = u*D*pBez;
end

% finds the x and y error of each point in "points" at the corresponding
% place on the Bezier curve
% INPUT:
%   points - matrix of points to fit as row vectors
%   pBez - matrix of Bezier control points as row vectors
% OUTPUT:
%   e - matrix of errors at each point as row vectors
function e = findErrors(points, pBez)
    n = size(points,1);
    u = 0:1/(n-1):1;
    e = zeros(size(points));
    
    for ii = 1:size(u,2)
        [p,~] = bezier(pBez,u(ii));
        e(ii,:) = p-points(ii,:);
    end
end

% finds the sum of the squared errors given a matrix of errors
% INPUT:
%   errors - matrix of positive and negative errors
% OUTPUT:
%   e - scalar value of sum of squared error
function e = sumSquaredError(errors)
    squaredErrors = errors.^2;
    e = sum(sum(squaredErrors));
end