xa = [175, 410, 675];
ya = [950, 2400, 1730];
xb = [160, 381, 656];
yb = [1008, 2500, 1760];
la = [60, 75, 42];
lb = [45, 88, 57];

xrot = [];
yrot = [];
fel = [];

for n = 1:3

    [x,y] = cir(xa(n), ya(n), la(n));
    plot(x,y)
    %axis("equal");
    hold on
    [x,y] = cir(xb(n), yb(n), lb(n));
    plot(x,y);

end

hold off
tol = 1e-12;

xstart = [200; 460; 700];
ystart = [1000; 2500; 1760];

for i = 1:3
    
    Xa = xa(i);
    Xb = xb(i);
    Ya = ya(i);
    Yb = yb(i);
    La = la(i);
    Lb = lb(i);
    Xstart = xstart(i);
    Ystart = ystart(i); 

    f = @(x,y) ((x-Xa).^2 + (y-Ya).^2 - La.^2);
    g = @(x,y) ((x-Xb).^2 + (y-Yb).^2 - Lb.^2);

    J =@(x,y) [2.*(x-Xa) 2.*(y-Ya); 2.*(x-Xb) 2.*(y-Yb)];

    dx = [1,1];

    while norm(dx) > tol
        A = [f(Xstart, Ystart); g(Xstart,Ystart)];
        B = [J(Xstart, Ystart)];

        h = -B\A;
        fel = [fel, h];
        dx = h;

        Xstart = Xstart + h(1);
        Ystart = Ystart + h(2);
    end

    xrot = [xrot, Xstart];
    yrot = [yrot, Ystart];

end

hold on

    plot(0,0,"o");
    plot(xrot(1), yrot(1), "o");
    plot(xrot(2), yrot(2), "o");
    plot(xrot(3), yrot(3), "o");
    plot(1020, 0, "o");

hold off

% fråga a
% --------------------
% Ai1 = [1 xrot(1) xrot(1).^2];
% Ai2 = [1 xrot(2) xrot(2).^2];
% Ai3 = [1 xrot(3) xrot(3).^2];
% 
% Ai = [Ai1; Ai2; Ai3];
% Bi = [yrot(1); yrot(2); yrot(3)];

%c = Ai\Bi;

% p = @(x) c(1) + c(2).*x + c(3).*x.^2;

% hold on
% 
% plot(u, p(u));
% 
% hold off

% fråga b
% ----------------------

Ai1 = [1 0 0 0 0];
Ai2 = [1 xrot(1) xrot(1).^2 xrot(1).^3 xrot(1).^4];
Ai3 = [1 xrot(2) xrot(2).^2 xrot(2).^3 xrot(2).^4];
Ai4 = [1 xrot(3) xrot(3).^2 xrot(3).^3 xrot(3).^4];
Ai5 = [1 1020 1020.^2 1020.^3 1020.^4];

Al = [Ai1; Ai2; Ai3; Ai4; Ai5];
Bl = [0; yrot(1); yrot(2); yrot(3); 0];

k = Al\Bl;

p = @(x) k(1) + k(2).*x + k(3).*x.^2 + k(4).*x.^3 + k(5).*x.^4;


u = 0:1:1020; 

hold on

plot(u, p(u));

hold off

%xdiff = abs(fel(2:end)-fel(1:end-1));

% normfel = sqrt(sum(fel.^2));
% disp(normfel)
% kon = abs(normfel(2:end)-normfel(1:end-1)); 

% längd av väg
% ------------------------

a = 0;
b = 1020;

n_int = 1000;

t = linspace(a,b,n_int+1);
step_length = (b-a)/n_int;

L = 0;

for j = 1:n
    L = L + (p(x(j)) + p(x(j+1)))/2*step_length;
end

disp(L)


