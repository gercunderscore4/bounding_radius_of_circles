function [center, radius] = get_circle(points)
    %[x, y, radius] = welzl(points, []);
    %center = [x, y];
    [radius, center] = ApproxMinBoundSphereND(points);
end