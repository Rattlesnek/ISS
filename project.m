% ISS - projekt
% Adam Pankuch (xpanku00)


pkg load signal


% SOURCE: https://www.mathworks.com/matlabcentral/answers/3058-plotting-circles
function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp);
end


% 1)

[s, Fs] = audioread('xpanku00.wav');
Fs                     % vzorkovacia frekvencia
N = length(s)          % dlzka vo vzorkoch
time = length(s) / Fs  % dlzka v sekundach
symb = length(s) / 16  % pocet sybolov (1 symbol 16 vzorkov)


% 2)

decode = zeros(symb, 1);

% kazdy 8 vzorok a podla neho dekodovane
j = 1;
for i = 8:16:(N - 8)
    if s(i) > 0 
        decode(j) = 1;
    end
    j++;
end

% TEST
%correct_symb = textread('xpanku00.txt', '%d');
%err_test = xor(correct_symb, decode);
%sum(err_test)

samples_s = 0.02 * Fs;
step = 0.02 / samples_s;
x_s = 0:step:(0.02 - step);

hold on;
% vykreslit prvych 20ms nacitaneho signalu
plot(x_s, s(1:samples_s));

samples_dec = samples_s / 16;
x_dec = (step*8):(step*16):(0.02 - step*8);
% vykreslit prvych 20ms dekodovaneho signalu
stem(x_dec, decode(1:samples_dec));
set(gca, 'Xtick', 0: 0.002 : 0.02);
set(gca, 'Ytick', -1: 0.2 : 1);

title('Uloha 2 - dekodovany signal');
xlabel('t[s]');
ylabel('s[n], symbols');
print -depsc decoded_sig.eps
hold off;


% 3)

% filter - prenosova funkcia
B = [0.0192 -0.0185 -0.0185 0.0192];
A = [1.0000 -2.8870 2.7997 -0.9113];

zplane(B, A);
title('Uloha 3 - filter');
xlabel('Real');
ylabel('Imag');
print -depsc filter.eps

nom = roots(B);
denom = roots(A);

hold on;
% nulove body
plot(real(nom), imag(nom), 'bo');
% poly
plot(real(denom), imag(denom), 'b*');
circle(0, 0, 1);
title('Uloha 3 - filter');
xlabel('Real');
ylabel('Imag');
print -depsc filter2.eps
hold off;


% 4)

% kmitoctova chcarakteristika filtru
H = freqz(B, A, Fs);
f = (0:length(H) - 1) / length(H) * Fs / 2;
plot(f, abs(H));
title('Uloha 4 - kimitoctova charakteristika filtru');
xlabel('f[Hz]');
ylabel('|H(z)|');
print -depsc filter_char.eps


% 5) 6)

shift = 17; % posunute o tuto hodnotu

ss = filter(B, A, s);
ss_shif = ss(shift:length(ss));

% dekodovanie filtrovaneho posunuteho signalu
filter_decode = zeros(length(ss_shif) / 16, 1);
j = 1;
for i = 8:16:(length(ss_shif) - 8)
    if ss_shif(i) > 0
        filter_decode(j) = 1;
    end
    j++;
end
    
hold on;
plot(x_s, s(1:samples_s));
% filtrovany signal
plot(x_s, ss(1:samples_s));
% filtrovany posunuty signal
plot(x_s, ss_shif(1:samples_s));
% dekodovany filtrovany posunuty signal
stem(x_dec, filter_decode(1:samples_dec));
set(gca, 'Xtick', 0: 0.002 : 0.02);
set(gca, 'Ytick', -1: 0.2 : 1);

title('Uloha 6 - filtrovany a posunuty signal');
xlabel('t[s]');
ylabel('s[n], ss[n], ss_{shifted}[n], symbols');
print -depsc filter_decode.eps
hold off; 


% 7)

err_vec = xor(filter_decode, decode(1:length(filter_decode)));
err_num = sum(err_vec) % pocet chyb

err_rate = err_num / length(filter_decode) % chybovost v %

% TEST
%err_test = xor(filter_decode(1:samples_dec), decode(1:samples_dec));
%sum(err_test)


% 8)

% spektrum signalu s
s_spec = fft(s);              
L = length(s); % dlzka signalu

s_spec_abs = abs(s_spec/L);
s_plot = s_spec_abs(1:L/2+1);
s_plot(2:end-1) = 2*s_plot(2:end-1);

f = Fs*(0:(L/2))/L;
hold on;
plot(f,s_plot); 

% spektrum signalu ss
ss_spec = fft(ss);
L = length(ss);

ss_spec_abs = abs(ss_spec/L);
ss_plot = ss_spec_abs(1:L/2+1);
ss_plot(2:end-1) = 2*ss_plot(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,ss_plot); 

title('Uloha 8 - Spektrum s a ss');
xlabel('f[Hz]');
ylabel('|S[f]|, |SS[f]|');
print -depsc spectrum.eps
hold off;


% 9)

% funkcia hustoty rozdelenia pravdepodobnosti
[s_pdf, x_pdf] = hist(s, 100);
s_norm = s_pdf / trapz(x_pdf, s_pdf);
plot(x_pdf, s_norm); % pdf
title('Uloha 9 - funkcia hustoty rozdelenia pravdepodobnosti');
print -depsc pdf.eps

integ_s = trapz(x_pdf, s_norm) % integral pdf = 1


% 10) 11)

% korelacne koeficienty
[R, lags] = xcorr(s, 50, 'biased');

plot(-50:50, R);
title('Uloha 10 - korelacne koeficienty (-50 ... 50)');
print -depsc correl.eps

R(51) % R[0]
R(52) % R[1]
R(67) % R[16]


% 12)

% zdruzena fce hustoty rozdelenia pravdepodobnosti
xmin = min(min(s)); 
xmax = max(max(s));
x = linspace(xmin,xmax,200);

s_plusOne = [s(2:end); 0];
[h,p,r] = hist2opt(s, s_plusOne, x);
imagesc (x,x,p);
axis xy;
colorbar;
title('Uloha 12 - zdruzena f. hustoty. roz. pravdep.');
xlabel('s[n+1]');
ylabel('s[n]');
print -depsc pdf2D.eps


% 13)

% 2 rozne sposoby vypoctu integralu
integ = trapz(x, trapz(x, p)) % 2 integral zdruz. funkc. pravdepodobnosti
surf = (x(2) - x(1))^2;
integ = sum(sum(p)) * surf % 2 integral zdruz. funkc. pravdepodobnosti


% 14)

% korelacny koeficient R[1] ziskany z funkcie hist2opt
r % R[1]





























