function [x, y, r] = three_point_circle(P)
    % solve circle equation for all points
    % x^2 + y^2 + A*x + B*y + C = 0
    % A*x + B*y + C*1 = -1*x^2 + -1*y^2
    % [x1, y1, 1] [A]   [(-x1*x1 -y1*y1)]
    % [x2, y2, 1] [B] = [(-x2*x2 -y2*y2)]
    % [x3, y3, 1] [C]   [(-x3*x3 -y3*y3)]
    % Ts = Q
    % s = T\Q
    %
    % convert equation to useful form
    % x^2 + y^2 + A*x + B*y = -C
    % (x^2 + A*x + (A/2)^2) + (y^2 + B*y + (B/2)^2) = (A/2)^2 + (B/2)^2 - C
    % (x + A/2)^2 + (y + B/2)^2 = (A/2)^2 + (B/2)^2 - C
    % (x - x0)^2 + (y - y0)^2 = r^2
    % x0 = -A/2
    % y0 = -B/2
    % r = (A/2)^2 + (B/2)^2 - C
    
    T = [P ones(3,1)];
    Q = -1*sum(P.^2, 2);
    if rank(T) < 3
        disp('GDI')
        P
        pause
        [x,y,r] = welzl(P,[]);
    else
        s = T\Q;
        
        x = -1*s(1)/2;
        y = -1*s(2)/2;
        r = sqrt(x^2 + y^2 - s(3));
    end
end
