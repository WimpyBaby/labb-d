function  f = felfunk(la,lb)

xa = [175, 410, 675];
ya = [950, 2400, 1730];
xb = [160, 381, 656];
yb = [1008, 2500, 1760];

x_guess = [];
y_guess = [];
xrot = [];
yrot = [];
fel = [];
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

Ai1 = [1 0 0 0 0];
Ai2 = [1 xrot(1) xrot(1).^2 xrot(1).^3 xrot(1).^4];
Ai3 = [1 xrot(2) xrot(2).^2 xrot(2).^3 xrot(2).^4];
Ai4 = [1 xrot(3) xrot(3).^2 xrot(3).^3 xrot(3).^4];
Ai5 = [1 1020 1020.^2 1020.^3 1020.^4];

Al = [Ai1; Ai2; Ai3; Ai4; Ai5];
Bl = [0; yrot(1); yrot(2); yrot(3); 0];

k = Al\Bl;

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
    f = TSim;
end

% nog = log2(abs(intvalues(2)-intvalues(1))/abs(intvalues(3)-intvalues(2))); %beräknar noggranhetsordningen
% 
% disp("Längden för vägen är : " + round(Trapets/1000, 2) + " kilometer")
% disp("Trapetsregeln har noggranhetsordning: " + round(nog))
end
	

