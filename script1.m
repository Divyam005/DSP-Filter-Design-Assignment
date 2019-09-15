m = 40;
q = ceil(0.1*m - 1);
r = m - 10*q;
bl = 5 + 1.4*q + 4*r;
bh = bl + 10;
fs = 320;

%% Normalization
wl = bl*2*pi/320;
wh = bh*2*pi/320;
s1 = wl - 2*2*pi/320;
s2 = wh + 2*2*pi/320;

%% Bilinear transform
wp1 = tan(wl/2);
wp2 = tan(wh/2);
ws1 = tan(s1/2);
ws2 = tan(s2/2);

%% Bandpass frequency transform
wo = sqrt(wp1*wp2);
B = wp2 - wp1;

wls1 = (ws1^2 - wo^2)/(B*ws1);
wls2 = (ws2^2 - wo^2)/(B*ws2);

%% Butterworth Filter general low pass
d1 = (1/(1-0.15)^2 - 1);
d2 = (1/0.15^2 - 1);
ws = min(abs(wls1),abs(wls2));
wp = 1;
N = ceil(log(d2/d1)/(2*log(ws/wp)));
wp/(d1^(1/(2*N)));
ws/(d2^(1/(2*N)));

%let wc be 1.08
wc = (ws/(d2^(1/(2*N)))+wp/(d1^(1/(2*N))))/2;
syms s z;
butter_eq = vpa(1+(s/(i*wc))^(2*N),4);
poles = solve(butter_eq,s);

%plot(real(poles),imag(poles),'x');

% Filter
poly = 1;
for j=1:N
    poly = poly*(s-poles(j));
end
H_analog(s) = (wc^N/expand(poly));
H_analog_bp(s) = H_analog((s^2+wo^2)/(B*s));
H_discrete(z) = H_analog_bp((z-1)/(z+1));

[nz, dz] = numden(vpa(H_discrete(z),5));                 
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          
k = dz(1);                                          
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 320e3);
plot(f,abs(H))
grid



