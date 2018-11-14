clear
clc
tic

%radii = [3 5 6 4*ones(1,15) 2*ones(1,5)]';
%radii = [3 4 5]';
%radii = [3, 4, 5]';
%radii = [3, 5, 6, 4*ones(1,15), 2*ones(1,5)]';
radii = rand(50,1);
%radii = ones(10,1);
radii = sort(radii, 'descend');
N = length(radii);

% calculate area in advance
areas = 4 * pi * (radii.^2);

% minimum separation (to avoid mathematical edge-cases)
e = 1E-3;

% calculate tightly packed positions
positions = zeros(N,2);
n = 32;
theta = linspace(0,2*pi,n)';
circ = [cos(theta), sin(theta)];
points = [];
if N > 1
    % circle 1 at origin
    % circle 2 beside it along x-axis
    positions(2,:) = [radii(1) + radii(2) + e, 0];

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
        d_ii = 0;
        d_jj = 0;
        d_mm = 0;
        d_collided = 0;
        d_a = 0;
        d_d = 0;
        d_h = 0;
        
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
                        collided = sum(sqrt(sum((positions(1:(nn-1),:) - positions(nn,:)).^2, 2)) < (radii(1:(nn-1)) + radii(nn)));
                        if collided
                            continue
                        end

                        % calculate new enclosing circle
                        new_points = [points; (positions(nn,:) + (radii(nn)*circ))];
                        [center, radius] = get_circle(new_points);
                        %center = (areas(1:nn)' * positions(1:nn,:)) / sum(areas(1:nn));
                        %radius = max(sqrt(sum((positions - center).^2, 2)) + radii);
    
                        % check if better
                        distance = norm(positions(nn,:) - center);
                        best = (radius < total_radius) || ((radius == total_radius) && (distance < dist_to_center));
                        
                        if best
                            total_radius = radius;
                            dist_to_center = distance;
                            pnp = positions(nn,:);
                            keep_points = positions(nn,:);
                            
                            % get/print debug for sub-steps
                            %nn;
                            %d_collided = collided;
                            %d_ii = ii;
                            %d_jj = jj;
                            %d_mm = mm;
                            %d_a = a;
                            %d_d = d;
                            %d_h = h;
                            %disp('')
                        end
            
                    end
                end
            end
        end
    
        % print debug for steps
        %nn
        %pnp
        %total_radius
        %dist_to_center
        %d_collided
        %d_ii
        %d_jj
        %d_mm
        %d_a
        %d_d
        %d_h
        %disp('')
        
        positions(nn,:) = pnp;
        points = [points; (positions(nn,:) + (radii(nn)*circ))];
    end
end

% re-calculate center and radius
[x,y,radius] = welzl(positions);
center = [x,y];
%center = (areas' * positions) / sum(areas);
%radius = max(sqrt(sum((positions - center).^2, 2)) + radii);

% print values
center
radius
N
usage = sum(areas) / (4 * pi * (radius^2))

figure(1);
clf;
plot_circles(positions, radii);
plot_circles(center, radius);
hold on
%for nn = 1:N
%    text(positions(nn,1), positions(nn,2), int2str(nn)) % print order
%    %text(positions(nn,1), positions(nn,2), sprintf('%.2f', radii(nn))) % print size
%end
%text(center(1), center(2), 'C')
plot(center(1), center(2), '*')
hold off
axis('equal');
toc
