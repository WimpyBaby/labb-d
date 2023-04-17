function [x,y] = cir(g,a,r)
t=linspace(0,2*pi);
x = g + r*cos(t);
y = a + r*sin(t);
end