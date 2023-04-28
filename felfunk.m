function  f = felfunk(la,lb)
xa = [175, 410, 675];
ya = [950, 2400, 1730];
xb = [160, 381, 656];
yb = [1008, 2500, 1760];
x_varden = [];
y_varden = [];
xrotter = [];
yrotter = [];
fel = [];
tol = 1e-12;
xstart = [200; 460; 700];
ystart = [1000; 2500; 1760];
for i = 1:3
    
    Xstart = xstart(i);
    Ystart = ystart(i); 
    f = @(x,y) ((x-xa(i)).^2 + (y-ya(i)).^2 - la(i).^2); %skriver om La funktionen till = 0 (funktionsberäkning)
    g = @(x,y) ((x-xb(i)).^2 + (y-yb(i)).^2 - lb(i).^2); % skriver om Lb funktionen till = 0 (funktionsberäkning)
    J =@(x,y) [2.*(x-xa(i)) 2.*(y-ya(i)); 2.*(x-xb(i)) 2.*(y-yb(i))]; %Jacobianen
    dx = [1,1];
    while norm(dx) > tol
        A = [f(Xstart, Ystart); g(Xstart,Ystart)];
        B = [J(Xstart, Ystart)];
        h = -B\A;
        fel = [fel, h];
        dx = h;
        Xstart = Xstart + h(1); % uppdatera x-startgissningen
        Ystart = Ystart + h(2); % uppdatera y-startgissningen
        x_varden = [x_varden, Xstart];
        y_varden = [y_varden, Ystart];
    end
    xrotter = [xrotter, Xstart];
    yrotter = [yrotter, Ystart];
end
%konvergens
% --------------------
% fråga b
% ----------------------
Ai1 = [1 0 0 0 0];
Ai2 = [1 xrotter(1) xrotter(1).^2 xrotter(1).^3 xrotter(1).^4];
Ai3 = [1 xrotter(2) xrotter(2).^2 xrotter(2).^3 xrotter(2).^4];
Ai4 = [1 xrotter(3) xrotter(3).^2 xrotter(3).^3 xrotter(3).^4];
Ai5 = [1 1020 1020.^2 1020.^3 1020.^4];
Al = [Ai1; Ai2; Ai3; Ai4; Ai5]; % A matrisen i Ay = b
Bl = [0; yrotter(1); yrotter(2); yrotter(3); 0]; %b matrisen ekvationen (alla funktionsvärden)
C = Al\Bl; %koefficientmatricen med C1 C2 C3 C4 C5 (y matrisen)
C1 = C(1);
C2 = C(2);
C3 = C(3);
C4 = C(4);
C5 = C(5);
p = @(x) C(1) + C(2).*x + C(3).*x.^2 + C(4).*x.^3 + C(5).*x.^4;
u = 0:1:1020; 
% längd av väg
% ------------------------
a = 0;
b = 1020;
n_int = [100 200 400];
p_length = @(x) sqrt(1 + (C(2) + 2.*C(3).*x + 3.*C(4).*x.^2 + 4.*C(5).*x.^3).^2); %formeln för kurvlängd: integral(sqrt(1 + f'(x)^2))
intvalues = [];
for w = n_int
    h1 = (b-a)/w;
    t1 = a:h1:b;
   Trapets = (h1/2)*(2*(sum(p_length(t1(2:end-1)))) + p_length(t1(1)) + p_length(t1(end)));
    intvalues = [intvalues, Trapets];
    f = Trapets;
end
% feltrapets = [(intvalues(1)-intvalues(2)); intvalues(2)-intvalues(3)]
% nog = log2(abs(intvalues(2)-intvalues(1))/abs(intvalues(3)-intvalues(2))); %beräknar noggranhetsordningen
% 
% disp("Längden för vägen är : " + round(Trapets/1000, 2) + " kilometer")
% disp("Trapetsregeln har noggranhetsordning: " + round(nog))
end
	

