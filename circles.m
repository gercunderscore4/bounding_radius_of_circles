clear
clc
tic

%%%%%%%%%%%%%%%% BEGIN INPUT %%%%%%%%%%%%%%%%

% NOTE: Diameter can be used in place of radius. Only impact is on the plot.

%radii = [3 5 6 4*ones(1,15) 2*ones(1,5)]';
%radii = [3 4 5]';
%radii = [3, 4, 5]';
%radii = [3, 5, 6, 4*ones(1,15), 2*ones(1,5)]';
radii = [0.142*ones(1,34)]';
%radii = ones(2,1);

%%%%%%%%%%%%%%%% END INPUT %%%%%%%%%%%%%%%%

% largest -> smallest gives good performance
radii = sort(radii, 'descend');
N = length(radii);

% minimum separation (to avoid mathematical edge-cases)
e = 1E-3;

% calculate tightly packed positions
positions = zeros(N,2);
% initialize values for convex hull
n = 8; % number of points along circle for convex hull calculation
theta = linspace(0,2*pi,n)';
circ = [cos(theta), sin(theta)];
points = bsxfun(@plus, positions(1,:), radii(1)*circ);
disp(sprintf('%d/%d', 1, N));
if N > 1
    % circle 1 at origin
    % circle 2 beside it along x-axis
    positions(2,:) = [radii(1) + radii(2) + e, 0];
    % maintain convex hull for fast enclosing circle calculations
    points = [points; bsxfun(@plus, positions(2,:), radii(2)*circ)];
    convex = convhull(points(:,1), points(:,2));
    points = points(convex(2:end),:); 
    disp(sprintf('%d/%d', 2, N));

    % place other circles touching two other circles
    %
    % x1^2 + y1^2 == (rn + r1 + e)^2
    % x2^2 + y2^2 == (rn + r2 + e)^2
    %
    % r1n = r1 + rn + e
    % r2n = r2 + rn + e
    %
    % https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
    %                           ___
    %                       ___/ | \_
    %                   ___/     |   \_
    %           r1  ___/         |     \_  r2
    %           ___/           h |       \_
    %       ___/                 |         \_
    %   ___/                    _|_          \_
    %  /_______________________|_|_|___________\
    %             a                     b
    %
    % r1 = r1n
    % r2 = r2n
    % d = sqrt((x1 - x2)^2 + (y1 - y2)^2) == a + b
    %
    % a = (r1^2 - r2^2 + d^2) / (2 * d)
    % h = sqrt(r1^2 - a^2)
    %
    % U = [x2 - x1, y2 - y1] / d
    % UN1 = [U(2), -U(1)]
    % UN2 = [-U(2), U(1)]
    %
    % xn = x1 + a*U(1) + h*UN(1)
    % yn = y1 + a*U(2) + h*UN(2)
    
    for nn = 3:N
        % get new radius
        rn = radii(nn);
    
        % clear optimization arrays
        total_radius   = Inf;
        dist_to_center = Inf;
        pnp = [NaN, NaN];
        
        for ii = 2:(nn-1)
            for jj = 1:(ii-1)
                % get two circles
                x1 = positions(ii,1);
                y1 = positions(ii,2);
                r1 = radii(ii);
        
                x2 = positions(jj,1);
                y2 = positions(jj,2);
                r2 = radii(jj);
        
                % solve triangle
                r1n = r1 + rn + e;
                r2n = r2 + rn + e;
                d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
                a = (r1n^2 - r2n^2 + d^2) / (2 * d);
                h = sqrt(r1n^2 - a^2);
    
                if (a >= 0) && (d >= 0) && (h >= 0) && (imag(a) == 0) && (imag(d) == 0) && (imag(h) == 0)
            
                    % unit vectors
                    U = [x2 - x1, y2 - y1] / d;
                    UN1 = [U(2), -U(1)];
                    UN2 = [-U(2), U(1)];
            
                    % potential new points
                    pn(1,:) = [x1 + a*U(1) + h*UN1(1), y1 + a*U(2) + h*UN1(2)];
                    pn(2,:) = [x1 + a*U(1) + h*UN2(1), y1 + a*U(2) + h*UN2(2)];
                    
                    for mm = 1:2
                        % get possible point
                        positions(nn,:) = pn(mm,:);
                        
                        % determine collision
                        collided = sum(sqrt(sum(bsxfun(@minus, positions(1:(nn-1),:), positions(nn,:)).^2, 2)) < (radii(1:(nn-1)) + radii(nn)));
                        if collided
                            continue
                        end

                        % maintain convex hull for fast enclosing circle calculations
                        new_points = [points; bsxfun(@plus, positions(nn,:), (radii(nn)*circ))];
                        convex = convhull(new_points(:,1), new_points(:,2));
                        new_points = new_points(convex(2:end),:); 
                        % calculate new enclosing circle
                        [center, radius] = get_circle(new_points);
    
                        % check if better
                        distance = norm(positions(nn,:) - center);
                        best = (radius < total_radius) || ((radius == total_radius) && (distance < dist_to_center));
                        
                        if best
                            total_radius = radius;
                            dist_to_center = distance;
                            pnp = positions(nn,:);
                            keep_points = positions(nn,:);
                        end
            
                    end
                end
            end
        end
        positions(nn,:) = pnp;
        % maintain convex hull for fast enclosing circle calculations
        points = [points; bsxfun(@plus, positions(nn,:), radii(nn)*circ)];
        convex = convhull(points(:,1), points(:,2));
        points = points(convex(2:end),:);
        disp(sprintf('%d/%d', nn, N));
    end
end

% re-calculate center and radius
[center,radius] = get_circle(points);

% re-center circles
points    = bsxfun(@minus, points   , center);
positions = bsxfun(@minus, positions, center);
center = [0,0];

% print values
clc
radius
N
usage = sum(radii.^2) / radius^2

figure(1);
clf;
plot_circles([positions; center], [radii; radius]);
axis('equal');
toc
