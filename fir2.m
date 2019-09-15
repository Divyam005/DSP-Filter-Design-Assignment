wp1 = 0.794194622827500;
ws1 = 0.844460105284936;
ws2 = 0.995256552657246;
wp2 = 1.045522035114683;
fs = 250e3;

%% To bring the passband inside 0.85 tolerance
wp1 = wp1*4/7 + ws1*3/7;
wp2 = wp2*4/7 + ws2*3/7;


%% Parameters for Kaiser window
delta_w = 2*pi*2e3/fs;
alpha = -20*log10(0.15);
if alpha<21
    beta = 0;
end
n_min = (alpha-7.95)/(2.285*delta_w);

N = ceil(n_min)+30;
r = rem(N,2);
if r == 1
else N = N + 1;
end

kaiser_win = kaiser(N,beta);

%% Ideal Bandstop filter in time
ideal_filter_h = zeros(N,1);
ideal_filter_l = zeros(N,1);
unity_filter = zeros(N,1);
p = (N+1)/2;
for i=1:N
    if i< (N+1)/2
        ideal_filter_h(i) = sin(wp2*(p-i))/(pi*(p-i));
        ideal_filter_l(i) = sin(wp1*(p-i))/(pi*(p-i));
        unity_filter(i) = sin(pi*(p-i))/(pi*(p-1));
    elseif i == (N+1)/2
        ideal_filter_h(i) = wp2/pi;
        ideal_filter_l(i) = wp1/pi;
        unity_filter(i) = 1;
    else
        ideal_filter_h(i) = sin(wp2*(i-p))/(pi*(i-p));
        ideal_filter_l(i) = sin(wp1*(i-p))/(pi*(i-p));
        unity_filter(i) = sin(pi*(i-p))/(pi*(i-p));
    end
end
fir_filter = (unity_filter - ideal_filter_h + ideal_filter_l).*kaiser_win;
fvtool(fir_filter);
[H,f] = freqz(fir_filter,1,1024, fs);
plot(f,abs(H));
grid
