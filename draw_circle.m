function draw_circle(x,y,r,P)
    figure(1)
    clf
    hold on
    % plot points
    if size(P,1) > 0
        plot(P(:,1), P(:,2), '*')
    end
    
    % plot circle
    if isnan(x)
    elseif r == 0
        plot(x, y, 'x')
    else
        t = linspace(0, 2*pi, 36);
        plot(x+r*cos(t), y+r*sin(t))
    end
    % keep it circular
    axis('equal')
    hold off
    
    disp('paused')
    pause
end
