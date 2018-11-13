function [x, y, r, P, R] = welzl(P, R)
    % https://en.wikipedia.org/wiki/Smallest-circle_problem#Welzl's_algorithm
    disp('NOT WORKING!')
    %disp('Welzl')
    %P
    %R
    %pause
    if (size(P,1) == 0) || (size(R,1) >= 3)
        %disp('B1')
        if size(R,1) == 1
            disp('B11')
            % P is empty, smallest "circle" containing the one point of R has radius zero
            p = R;
            r = 0;
            x = R(1,1);
            y = R(1,2);
            r = 0;
        elseif size(R,1) == 2
            disp('B12')
            % P is empty, R has two points, smallest circle is centered at midpoint
            p = [R(1,:) + R(2,:)] / 2;
            x = p(1);
            y = p(2);
            r = norm(R(1,:) - R(2,:), 2) / 2;
        elseif size(R,1) == 3
            disp('B13')
            % return the circle defined by the tree points
            % x^2 + y^2 + A*x + B*y + C = 0
            % A*x + B*y + C*1 = -1*x^2 + -1*y^2
            % [x1, y1, 1] [A]   [(-x1*x1 -y1*y1)]
            % [x2, y2, 1] [B] = [(-x2*x2 -y2*y2)]
            % [x3, y3, 1] [C]   [(-x3*x3 -y3*y3)]
            % Ps = Q
            % s = P\Q
            %
            % x^2 + y^2 + A*x + B*y = -C
            % (x^2 + A*x + (A/2)^2) + (y^2 + B*y + (B/2)^2) = (A/2)^2 + (B/2)^2 - C
            % (x + A/2)^2 + (y + B/2)^2 = (A/2)^2 + (B/2)^2 - C
            % x = -A/2
            % y = -B/2
            % r = sqrt((A/2)^2 + (B/2)^2 - C)
            
            P = [R ones(3,1)];
            Q = -1*sum(R.^2, 2);
            s = P\Q;
            
            x = -1*s(1)/2;
            y = -1*s(2)/2;
            r = sqrt(x^2 + y^2 - s(3));
        else
            disp('B14')
            x = NaN;
            y = NaN;
            r = NaN;
        end
    else
        %disp('B2')
        i = randi(size(P,1));
        p = P(i,:);
        P(i,:) = []; % delete row
        [x, y, r, ~, ~] = welzl(P, R);
        
        if isnan(r) || (norm([p; [x, y]], 2) > r)
            %disp('B21')
            [x, y, r, ~, ~] = welzl(P, [R; p]);
        end
    end
    %x
    %y
    %r
end
