function [x,y] = boll(g,a,r)
t=linspace(0,2*pi);
x = g + r*cos(t);
y = a + r*sin(t);
end