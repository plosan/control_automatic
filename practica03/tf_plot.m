clear; close all; clc;

T0 = 1e-3;
K = 1;
tau = 1/10;

G = @(s) K./(tau*s + 1).*exp(-T0*s);

g = @(s) 20*log10(abs(G(1i*s)));

f = @(s) K./sqrt(1 + tau^2*s.^2);
f2 = @(s) K./(tau*s);

% Interpreter latex
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure(1);
hold on;
fplot(g, [1e-6 1e6], 'b');
yline(20*log10(K));

set(gca, 'XScale', 'log');
xlim([1e-2 1e2]);
ylim([0 1.2]);
xlabel('Freq\"u\`encia $\omega \ [\mathrm{rad}/\mathrm{s}]$');
ylabel('$20 \log_{10} |G(s)|$');
box on;

hold off;

wc = 1/tau;

h = figure(2);
hold on;
fplot(@(s) 20*log10(f(s)), [1e-6 1e6], 'b', 'LineWidth', 1);
fplot(@(s) 20*log10(f2(s)), [1e-6 1e6], 'r');
fplot(@(s) 20*log10(f2(s)), [wc 1e6], '--k');
plot([1e-6 wc], [20*log10(K) 20*log10(K)], '--k');
% plot([1e-6 wc], [20*log10(K)-20*log10(sqrt(2)) 20*log10(K)-20*log10(sqrt(2))], '--k');
% plot([wc wc], [20*log10(K) -50], '--k');
% scatter(wc, 20*log10(K)-20*log10(sqrt(2)), 30, 'r', 'filled');
% text(0.3, 20*log10(K)-20*log10(sqrt(2))+1, '$20 \log_{10} K - 20 \log_{10} \sqrt{2}$', 'HorizontalAlignment', 'center', 'FontSize', 11);
% text(wc+2.5, -10, '$\omega_c$', 'HorizontalAlignment', 'center', 'FontSize', 11);

text(0.3, 20*log10(K)+0.7, '$20 \log_{10} K$', 'HorizontalAlignment', 'center', 'FontSize', 11);

% yline(20*log10(K), '--k');
% xline(wc, '--k');
xlim([1e-2 1e2]);
ylim([-20 5]);
xlabel('Freq\"u\`encia $\omega$');
ylabel('M\`odul');
yticks([]);
xticklabels({});
ytickangle(90);
legend("$g(\omega) = 20 \log_{10}|G(i \omega)|$", "$h(\omega) = 20 \log_{10} (K / \tau) - 20 \log_{10} \omega$", ...
    "Location", "northwest");
set(gca, 'XScale', 'log');
set(gca, 'Xgrid', 'on');


set(gcf, 'units', 'centimeters', 'position', [22,5,17,10]);
hold off;

save2pdf(h, "tf_plot.pdf");



h = figure(3);
hold on;
title('\textbf{M\`odul -- Freq\"u\`encia}');
fplot(@(s) 20*log10(f(s)), [1e-6 1e6], 'b', 'LineWidth', 1);
fplot(@(s) 20*log10(f2(s)), [1e-6 1e6], 'r');
fplot(@(s) 20*log10(f2(s)), [wc 1e6], '--k');
plot([1e-6 wc], [20*log10(K) 20*log10(K)], '--k');
plot([1e-6 wc], [20*log10(K)-20*log10(sqrt(2)) 20*log10(K)-20*log10(sqrt(2))], '--k');
plot([wc wc], [20*log10(K) -50], '--k');
scatter(wc, 20*log10(K)-20*log10(sqrt(2)), 30, 'r', 'filled');

text(0.3, 20*log10(K)+1, '$20 \log_{10} K$', 'HorizontalAlignment', 'center', 'FontSize', 11);
text(0.3, 20*log10(K)-20*log10(sqrt(2))+1, '$20 \log_{10} K - 20 \log_{10} \sqrt{2}$', 'HorizontalAlignment', 'center', 'FontSize', 11);
text(wc+2.5, -10, '$\omega_c$', 'HorizontalAlignment', 'center', 'FontSize', 11);


% yline(20*log10(K), '--k');
% xline(wc, '--k');
xlim([1e-2 1e2]);
ylim([-20 5]);
xlabel('Freq\"u\`encia $\omega$');
ylabel('M\`odul [dB]');
yticks([]);
xticklabels({});
ytickangle(90);
legend("$g(\omega) = 20 \log_{10}|G(i \omega)|$", "$h(\omega) = 20 \log_{10} (K / \tau) - 20 \log_{10} \omega$", ...
    "Location", "southwest");
set(gca, 'XScale', 'log');
set(gca, 'Xgrid', 'on');

box on;
set(gcf, 'units', 'centimeters', 'position', [22,5,17,8]);
hold off;

save2pdf(h, "tf_plot_2.pdf");
