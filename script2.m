m = 40;
q = ceil(0.1*m - 1);
r = m - 10*q;
bl = (5 + 1.2*q + 2.5*r)*1e3;
bh = bl + 6e3;
fs = 250e3;

%% Normalization
wl = bl*2*pi/fs;
wh = bh*2*pi/fs;
p1 = wl - 2*2e3*pi/fs;
p2 = wh + 2*2e3*pi/fs;

%% Bilinear transform
ws1 = tan(wl/2);
ws2 = tan(wh/2);
wp1 = tan(p1/2);
wp2 = tan(p2/2);

%% Bandstop - Lowpass frequency transform
wo = sqrt(wp1*wp2);
B = wp2 - wp1;

wls1 = (B*ws1)/(ws1^2 - wo^2);
wls2 = (B*ws2)/(ws2^2 - wo^2);

%% Chebyschev Filter parameters for low pass
d1 = (1/(1-0.15)^2 - 1);
d2 = (1/0.15^2 - 1);

epsilon = sqrt(d1);

ws = min(abs(wls1),abs(wls2));
wp = 1;

n_min = acosh(sqrt(d2/d1))/acosh(sqrt(ws/wp));
N = ceil(n_min)+1;

syms s z;
cheb_eq = vpa(1 + (epsilon^2)*(cos(N*acos(s/(i*wp))))^2,4);
poles = (solve(cheb_eq,s));

%plot(real(poles),imag(poles),'x');

% Filter
poly = 1;
rpoles = 1;
for k = 1:2*N
    if real(poles(k)) < 0
        poly = poly*(s-poles(k));
        rpoles = rpoles*abs(poles(k));
    end
end
H_analog(s) = rpoles*(-1)^N/(sqrt(1+d1)*expand(poly));
H_analog_bs(s) = H_analog((B*s)/(s^2+wo^2));
H_discrete(z) = H_analog_bs((z-1)/(z+1));


[nz, dz] = numden(vpa(H_discrete(z),4));                 
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                         
k = dz(1);                                         
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                      

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, fs);
plot(f,abs(H))
grid

