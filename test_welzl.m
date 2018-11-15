clear
clc

% get points
P = [0 1; 0 2; 0 3]
% calculate smallest enclosing circle
% https://en.wikipedia.org/wiki/Smallest-circle_problem#Welzl's_algorithm
%[x,y,r] = welzl(P, []);
[x,y,r] = three_point_circle(P, []);

figure(1)
clf
hold on
% plot circle
t = linspace(0, 2*pi, 36);
plot(x+r*cos(t), y+r*sin(t))
% plot points
plot(P(:,1), P(:,2), '*')
% keep it circular
axis('equal')
hold off
