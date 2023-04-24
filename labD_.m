xa = [175, 410, 675];
ya = [950, 2400, 1730];
xb = [160, 381, 656];
yb = [1008, 2500, 1760];
la = [60, 75, 42];
lb = [45, 88, 57];

x_guess = [];
y_guess = [];
xrot = [];
yrot = [];
fel = [];

for n = 1:3

    [x,y] = boll(xa(n), ya(n), la(n));
    plot(x,y)
    %axis("equal");
    hold on
    [x,y] = boll(xb(n), yb(n), lb(n));
    plot(x,y);

end

hold off
tol = 1e-12;

xstart = [200; 460; 700];
ystart = [1000; 2500; 1760];

for i = 1:3
    
    Xstart = xstart(i);
    Ystart = ystart(i); 

    f = @(x,y) ((x-xa(i)).^2 + (y-ya(i)).^2 - la(i).^2);
    g = @(x,y) ((x-xb(i)).^2 + (y-yb(i)).^2 - lb(i).^2);

    J =@(x,y) [2.*(x-xa(i)) 2.*(y-ya(i)); 2.*(x-xb(i)) 2.*(y-yb(i))];

    dx = [1,1];

    while norm(dx) > tol
        A = [f(Xstart, Ystart); g(Xstart,Ystart)];
        B = [J(Xstart, Ystart)];

        h = -B\A;
        fel = [fel, h];
        dx = h;

        Xstart = Xstart + h(1);
        Ystart = Ystart + h(2);

        x_guess = [x_guess, Xstart];
        y_guess = [y_guess, Ystart];
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

%konvergens
% --------------------

x_diff = abs(x_guess(2:end)-x_guess(1:end-1));
y_diff = abs(y_guess(2:end)-y_guess(1:end-1));

x_kon = abs(x_diff(2:end))./abs((x_diff(1:end-1))).^2;
y_kon = abs(y_diff(2:end))./abs((y_diff(1:end-1))).^2;

iter = [1:1:15]';

x_list = x_kon';
y_list = y_kon';

disp([iter, x_list, y_list])

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


% längd av väg
% ------------------------

a = 0;
b = 1020;
n_int = [100 200 400];

p_length = @(x) sqrt(1 + (k(2) + 2.*k(3).*x + 3.*k(4).*x.^2 + 4.*k(5).*x.^3).^2);

intvalues = [];

for w = n_int
    h1 = (b-a)/w;
    t1 = a:h1:b;
    TSim = (h1/2)*(2*(sum(p_length(t1(2:end-1)))) + p_length(t1(1)) + p_length(t1(end)));
    intvalues = [intvalues, TSim];
end

nog = log2(abs(intvalues(2)-intvalues(1))/abs(intvalues(3)-intvalues(2)));

disp("Längden för vägen är : " + round(TSim/1000, 2) + " kilometer")
disp("Trapetsregeln har noggranhetsordning: " + round(nog))


%-----------------------





