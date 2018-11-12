function plot_circles(centers, radii)
    t = linspace(0,2*pi,36);
    hold on;
    for ii = 1:length(radii)
        x = centers(ii,1);
        y = centers(ii,2);
        r = radii(ii);
        plot(x + (r * cos(t)), y + (r * sin(t)))
    end
    hold off;
end
