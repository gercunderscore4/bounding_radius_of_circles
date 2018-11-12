clear
clc

N = 10;
radii = [1 2 3 4 5 6 7 8 9 10]';
%radii = rand(N,1);
%radii = ones(N,1);
areas = 4 * pi * (radii.^2);
total_area = sum(areas);

% separation
e = 1E-3;

positions = zeros(N,2);

r1 = radii(1);
r2 = radii(2);
r3 = radii(3);
positions(1,:) = [0, 0];
positions(2,:) = [r1 + r2 + e, 0];

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

x1 = positions(1,1)
y1 = positions(1,2)
x2 = positions(2,1)
y2 = positions(2,2)
rn = r3
r1n = r1 + rn + e
r2n = r2 + rn + e
d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
a = (r1n^2 - r2n^2 + d^2) / (2 * d)
h = sqrt(r1n^2 - a^2)
U = [x2 - x1, y2 - y1] / d
UN1 = [U(2), -U(1)]
UN2 = [-U(2), U(1)]
xn = x1 + a*U(1) + h*UN1(1)
yn = y1 + a*U(2) + h*UN1(2)
positions(3,:) = [xn, yn];

for nn = 4:N
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
                    % try first
                    positions(nn,:) = pn(mm,:);
                    center = (areas(1:nn)' * positions(1:nn,:)) / sum(areas(1:nn))

                    % determine collision
                    collided = 0;
                    for kk = 1:(nn-1)
                        if (norm(positions(kk,:) - positions(nn,:))) < (radii(kk) + radii(nn))
                            collided = 1;
                            break
                        end
                    end

                    % calculate new total radius
                    radius = 0;
                    for kk = 1:nn
                        new_radius = norm(positions(kk,:) - center) + radii(kk);
                        if new_radius > radius
                            radius = new_radius;
                        end
                    end
                
                    % check if better
                    distance = norm(positions(nn,:) - center);

                    best = (radius < total_radius) || ((radius == total_radius) && (distance < dist_to_center));

                    %figure(1);
                    %clf;
                    %plot_circles(positions(1:nn,:), radii(1:nn));
                    %axis('equal');
                    %nn
                    %ii
                    %jj
                    %collided
                    %best
                    %pause;
        
                    if (collided == 0) && best
                        total_radius = radius;
                        dist_to_center = distance;
                        pnp = positions(nn,:);
                        nn;
                        d_collided = collided;
                        d_ii = ii;
                        d_jj = jj;
                        d_mm = mm;
                        d_a = a;
                        d_d = d;
                        d_h = h;
                        disp('')
                    end
        
                end
            end
        end
    end

    nn
    pnp
    total_radius
    dist_to_center
    d_collided
    d_ii
    d_jj
    d_mm
    d_a
    d_d
    d_h
    disp('')
    
    positions(nn,:) = pnp;
end

center = (areas(1:nn)' * positions(1:nn,:)) / sum(areas(1:nn))

% calculate new total radius
radius = 0;
for kk = 1:nn
    new_radius = norm(positions(kk,:) - center) + radii(kk);
    if new_radius > radius
        radius = new_radius;
    end
end


figure(1);
clf;
plot_circles(positions(1:nn,:), radii(1:nn));
plot_circles(center, radius);
hold on
for nn = 1:N
    text(positions(nn,1), positions(nn,2), int2str(nn))
end
text(center(1), center(2), 'C')
hold off
axis('equal');
