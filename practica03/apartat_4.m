clear; close all; clc;

%% 1. DADES
load("practica_03.mat");

plot1 = 0;
plot2 = 1;
plot31 = 0;
plot_modul_frequencia = 0;
plot_fase_frequencia = 0;
plot_tau = 0;
plot5 = 0;

% Interpreter latex
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 2. CALCULS
logmod20 = 20*log10(modul);

A1 = logmod20(1);   % 20*log(mod(G)) at the lowest frequency
K = 10^(A1/20);     % Proportionality constant
a = A1 - 3.01;      % To find cut frequency
j = 1;              % Find index of logmod20 that sandwiches a
while ((logmod20(j)-a)*(logmod20(j+1)-a)) > 0
    j = j + 1;
end
wc = (w(j+1) - w(j))*(a - logmod20(j))/(logmod20(j+1) - logmod20(j)) + w(j);    % Cut frequency interpolation
tau = 1/wc;         % Time parameter

T0 = -pi/180*(fase(end)+90)/w(end);

%% 3. MARGES TEORICS
G = tf([K], [tau 1], 'inputdelay', T0); % Theoretical model
[mag,phase,wout] = bode(G, w);          % Bode plot
mag = reshape(mag, [1 9]);              % Bode plot magnitude
phase = reshape(phase, [1 9]);          % Bode plot phase

Gf = @(omega) K/(1 + 1i*tau*omega).*exp(1i*T0*omega);
absGf = @(omega) abs(Gf(omega));
modGf = @(omega) 20*log10(absGf(omega));
angGf = @(omega) 180/pi*angle(Gf(omega));

h = figure();
margin(G)
set(gcf, 'units', 'centimeters', 'position', [22,5,17,12]);
exportgraphics(h, 'marges_teorics.pdf', 'ContentType', 'vector');

%% 4. MARGES EMPIRICS

% Extrapolació de mòdul i fase
mextra = @(omega) (logmod20(9)-logmod20(6))/(log(w(9))-log(w(6)))*(log(omega) - log(w(6))) + logmod20(6);
fextra = @(omega) (fase(9)-fase(5))/(log(w(9))-log(w(5)))*(log(omega)-log(w(5))) + fase(5);

% Freqüència on la fase és -180
wcg = fzero(@(omega) fextra(omega)+180, 1e3);

% Interpolació lineal per trobar quan el 20*log10(abs(G)) és 0
minter = @(omega) (logmod20(6)-logmod20(5))/(log(w(6))-log(w(5)))*(log(omega) - log(w(5))) + logmod20(5);
wcf = fzero(minter, 10);    % Freqüència on 20*log10(abs(G)) = 0

h = figure();
sgtitle("\textbf{Diagrama de Bode de la funci\'o de transfer\`encia emp\'irica}", "Fontsize", 11);
subplot(2,1,1);
hold on;    
plot(w, logmod20, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
fplot(mextra, [w(6) 1e4], 'r');
% Rectes per marge de guany
plot([1e-1 wcg wcg], [modGf(wcg) modGf(wcg) -200], '--k');
plot([wcg wcg], [modGf(wcg) 0], 'g', 'LineWidth', 2);
scatter(wcg, modGf(wcg), 30, 'g', 'filled');
% Rectes per marge de fase
plot([wcf wcf], [0 -200], '--k');
plot([1e-1 1e4], [0 0], '--k');
scatter(wcf, 0, 30, 'g', 'filled');
% Altres
xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
ylabel('M\`odul $20 \log_{10} |G(i \omega)| \ [\mathrm{dB}]$');
xlim([1e-1 1e4]);
set(gca, 'XScale', 'log');
grid on;
box on;
legend("M\`odul de la funci\'o de transfer\`encia emp\'irica", "Extrapolaci\'o lineal del m\`odul", "Location", "southwest");
hold off;

subplot(2,1,2);
hold on;
plot(w, fase, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
fplot(fextra, [w(5) 1e4], 'r');
% Rectes per marge de guany
plot([1e-1 wcg wcg], [-180 -180 180], '--k');
scatter(wcg, -180, 30, 'g', 'filled');
% Rectes per marge de fase
plot([wcf wcf 1e-1], [0 fextra(wcf) fextra(wcf)], '--k');
plot([wcf wcf], [fextra(wcf) 0], 'g', 'LineWidth', 2);
scatter(wcf, fextra(wcf), 30, 'g', 'filled');
% Altres
xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
ylabel('Fase $\theta \ [\mathrm{deg}]$');
xlim([1e-1 1e4]);
ylim([-200 0]);
set(gca, 'XScale', 'log');
yticks([-180:60:0]);
grid on;
box on;
hold off;
set(gcf, 'units', 'centimeters', 'position', [22,5,17,17]);
legend("Fase de la funci\'o de transfer\`encia emp\'irica", "Extrapolaci\'o lineal de la fase", "Location", "southwest");

save2pdf(h, "bode_empiric_gran.pdf");