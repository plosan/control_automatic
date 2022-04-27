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


%% 2. PLOTS ENTRADA I SORTIDA
% Interpreter latex
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

if plot1 == 1
    h1 = figure(1);
    sgtitle("\textbf{Entrada i sortida del sistema}", "Fontsize", 11);
    hold on;
    
    % Entrada
    subplot(3,2,1);
    plot(temps, entrada, 'b');
    xlabel("Temps $t \ [\mathrm{s}]$");
    ylabel("Multisinus entrada");
    grid on;
    box on;
    
    subplot(3,2,3);
    semilogx(w, gentrada, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('M\`odul entrada [dB]');
    ylim([0.40 0.50]);
    grid on;
    box on;
    
    subplot(3,2,5);
    semilogx(w, fentrada, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('Fase entrada [deg]');
    grid on;
    box on;
    
    % Sortida
    subplot(3,2,2);
    plot(temps, sortida, 'b');
    xlabel("Temps $t \ [\mathrm{s}]$");
    ylabel("Sortida");
    grid on;
    box on;
    
    subplot(3,2,4);
    semilogx(w, gsortida, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('M\`odul sortida [dB]');
    grid on;
    box on;
    
    subplot(3,2,6);
    semilogx(w, fsortida, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('Fase sortida [deg]');
    grid on;
    box on;
    
    set(gcf, 'units', 'centimeters', 'position', [5,5,17,15]);
    
    hold off;
    
    save2pdf(h1, 'entrada_sortida.pdf');
end

%% 3. DIAGRAMA DE BODE FUNCIO TRANSFERENCIA
logmod20 = 20*log10(modul);

if plot2 == 1

    h2 = figure(2);
    sgtitle("\textbf{Diagrama de Bode de la funci\'o de transfer\`encia emp\'irica}", "Fontsize", 11);
    hold on;    
    subplot(2,1,1);
    semilogx(w, logmod20, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('M\`odul $20 \log_{10} |G(i \omega)| \ [\mathrm{dB}]$');
    xlim([0.1 200]);
    ylim([-25 5]);
    grid on;
    box on;
    
    subplot(2,1,2);
    semilogx(w, fase, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('Fase $\theta \ [\mathrm{deg}]$');
    xlim([0.1 200]);
    ylim([-120 0]);
    grid on;
    box on;
    
    set(gcf, 'units', 'centimeters', 'position', [22,5,17,12]);
    hold off;
    
    save2pdf(h2, 'bode_funcio_transferencia.pdf');
    
    for i = 1:length(modul)
        fprintf("$%2d$ %3s $%7.3f$ %3s $%7.3f$ %3s $%7.3f$ \\\\ \n", i, "&", w(i), "&", 20*log10(modul(i)), "&", fase(i));
    end

end

%% 4. CALCULS

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

% Plot K constant
if plot_modul_frequencia == 1

    h_modul_frequencia = figure();
    hold on;
    title('\textbf{M\`odul -- Freq\"u\`encia}');
    plot(w, 20.*log10(modul), '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
%     scatter(wc, A1-3.01, 30, 'r', 'filled');
    set(gca, 'XScale', 'log');
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('M\`odul $20 \log_{10} | G(s) | \ [\mathrm{dB}]$');
    xlim([0.1 200]);
    ylim([-25 5]);
    grid on;
    box on;    
    legend("M\`odul de la funci\'o de transfer\`encia");
    set(gcf, 'units', 'centimeters', 'position', [5,5,17,8]);
    hold off;

end

% Plot pure delay
if plot_fase_frequencia == 1
    h_fase_frequencia = figure();
    title('\textbf{Fase -- Freq\"u\`encia}');
    hold on;
    plot(w, fase, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    set(gca, 'XScale', 'log');
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('Fase $\theta \ [\mathrm{deg}]$');
    xlim([0.1 200]);
    ylim([-120 0]);
    grid on;
    box on;
    legend("Fase de la funci\'o de transfer\`encia");
    set(gcf, 'units', 'centimeters', 'position', [22,5,17,8]);
    hold off;

end

if plot_tau == 1

    h_plot_tau = figure();
    hold on;
    title('\textbf{M\`odul -- Freq\"u\`encia}');
    plot(w, 20.*log10(modul), '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    
    yline(A1, '--k');
    yline(A1-3.01, '--k');
    xline(wc, '--k');
    
    text(60, A1+0.7, '$1.324 \ \mathrm{dB}$', 'HorizontalAlignment', 'center');
    text(60, A1-3.01+0.7, '$-1.686 \ \mathrm{dB}$', 'HorizontalAlignment', 'center');
    text(10+2, -17.5, '$\omega_c = 10.986 \ \mathrm{rad} / \mathrm{s}$');
    
    set(gca, 'XScale', 'log');
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('M\`odul $20 \log_{10} |G(i \omega)| \ [\mathrm{dB}]$');
    
    scatter(wc, A1-3.01, 30, 'r', 'filled');
    
    xlim([0.1 200]);
    ylim([-25 5]);
    grid on;
    box on;
    
    set(gcf, 'units', 'centimeters', 'position', [5,5,17,8]);
    
    legend("M\`odul de la funci\'o de transfer\`encia", "Location", "southwest");
    hold off;

    save2pdf(h_plot_tau, 'funcio_transferencia_parametres.pdf');

end


%% 5. SIMULAR FUNCIO TRANSFERENCIA
G = tf([K], [tau 1], 'inputdelay', T0); % Theoretical model
[mag,phase,wout] = bode(G, w);          % Bode plot
mag = reshape(mag, [1 9]);              % Bode plot magnitude
phase = reshape(phase, [1 9]);          % Bode plot phase

Gf = @(omega) K/(1 + 1i*tau*omega).*exp(1i*T0*omega);
absGf = @(omega) abs(Gf(omega));
angGf = @(omega) 180/pi*angle(Gf(omega));

f = @(omega) angGf(omega) + 180;

% Bode plot
if plot5 == 1

    h5 = figure(5);
    sgtitle("\textbf{Comparaci\'o dels diagrames de Bode emp\'iric i te\`oric}", "Fontsize", 11);
    subplot(2,1,1);
    hold on;
    plot(w, 20*log10(modul), '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    plot(wout, 20*log10(mag), 'r');
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('M\`odul $20 \log_{10} |G(i \omega)| \ [\mathrm{dB}]$');
    xlim([0.1 200]);
    ylim([-25 5]);
    grid on;
    box on;
    set(gca, 'XScale', 'log');
    legend("Model emp\'iric", "Model te\`oric");
    hold off;
    
    subplot(2,1,2);
    hold on;
    plot(w, fase, '-ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    plot(wout, phase, 'r');
    xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
    ylabel('Fase [deg]');
    xlim([0.1 200]);
    ylim([-120 0]);
    grid on;
    box on;
    set(gca, 'XScale', 'log');
    legend("Model emp\'iric", "Model te\`oric");
    hold off;
    
    set(gcf, 'units', 'centimeters', 'position', [22,5,17,12]);
    
    
    save2pdf(h5, 'bode_teoric_empiric.pdf');

end

%% 6. GAIN MARGIN AND PHASE MARGIN
% h6 = figure(6);
% margin(G);


[Gm, Pm, Wcg, Wcp] = margin(G);
fprintf("%10s = %7.3f dB at %.2f rad/s\n", "Gm", Gm, Wcg);
fprintf("%10s = %7.3f deg at %.2f rad/s\n", "Pm", Pm, Wcp);

figure();
margin(G);








