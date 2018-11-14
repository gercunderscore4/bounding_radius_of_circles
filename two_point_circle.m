function [x, y, r] = two_point_circle(P)
    % center at midpoint
    % diameter is distance between points
    p = [P(1,:) + P(2,:)] / 2;
    x = p(1);
    y = p(2);
    r = norm(P(1,:) - P(2,:), 2) / 2;
end
