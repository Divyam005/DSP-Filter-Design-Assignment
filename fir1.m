ws1 = 0.926769832808989;
wp1 = 0.966039740978862;
wp2 = 1.162389281828224;
ws2 = 1.201659189998096;
fs = 320e3;

%% To bring the passband inside 0.85 tolerance
wp1 = wp1*4/8 + ws1*4/8;
wp2 = wp2*4/8 + ws2*4/8;


%% Parameters for Kaiser window
delta_w = 2*pi*2e3/fs;
alpha = -20*log10(0.15);
if alpha<21
    beta = 0;
end
n_min = (alpha-7.95)/(2.285*delta_w);

N = ceil(n_min)+35;
r = rem(N,2);
if r == 1
else N = N + 1;
end

kaiser_win = kaiser(N,beta);

%% Time domain ideal filter
ideal_filter_h = zeros(N,1);
ideal_filter_l = zeros(N,1);
p = (N+1)/2;
for i=1:N
    if i< (N+1)/2
        ideal_filter_h(i) = sin(wp2*(p-i))/(pi*(p-i));
        ideal_filter_l(i) = sin(wp1*(p-i))/(pi*(p-i));
    elseif i == (N+1)/2
        ideal_filter_h(i) = wp2/pi;
        ideal_filter_l(i) = wp1/pi;
    else
        ideal_filter_h(i) = sin(wp2*(i-p))/(pi*(i-p));
        ideal_filter_l(i) = sin(wp1*(i-p))/(pi*(i-p));
    end
end
fir_filter = (ideal_filter_h-ideal_filter_l).*kaiser_win;

fvtool(fir_filter);
[H,f] = freqz(fir_filter,1,1024, fs);
plot(f,abs(H))
grid