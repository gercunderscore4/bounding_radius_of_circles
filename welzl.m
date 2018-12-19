function [x, y, r] = welzl(P, R)
    sP = size(P,1);
    sR = size(R,1);
    disp(sprintf('P = %4d    R = %4d', sP, sR))
    if sR > 3
        [x,y,r] = welzl(R, []);
    elseif (sP == 0)
        if sR == 3
            [x,y,r] = three_point_circle(R);
            %draw_circle(x,y,r,[P;R])
        elseif sR == 2
            [x,y,r] = two_point_circle(R);
            %draw_circle(x,y,r,[P;R])
        elseif sR == 1
            [x,y,r] = one_point_circle(R);
            %draw_circle(x,y,r,[P;R])
        else
            x = NaN;
            y = NaN;
            r = NaN;
            %draw_circle(x,y,r,[P;R])
        end
    else
        % take one point
        i = randi(sP);
        p = P(i,:);
        P(i,:) = [];
        % try building a circle without it
        [x, y, r] = welzl(P, R);
        
        % if point is in circle, it's not necessary, and this circle is good
        if point_in_circle(p, x, y, r);
            x = x;
            y = y;
            r = r;
        else
            % else P is necessary
            R = [R; p];
            [x,y,r] = welzl(P, R);
        end
    end
end
