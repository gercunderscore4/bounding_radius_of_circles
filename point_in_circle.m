function b = point_in_circle(p, x, y, r)
    if isnan(x) || isnan(y) || isnan(r)
        % circle not defined
        % can't be in an undefined circle
        b = 0;
    else
        % get distance
        d = norm((p - [x, y]), 2);
        % if distance less than radius
        b = d <= r;
end
