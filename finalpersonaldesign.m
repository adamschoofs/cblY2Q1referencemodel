%% ===================== PART C — Single plot (seconds) ===================
figure('Color','w'); hold on; grid on; box on;

% Model: tank water temperature (°C)
plot(t_sim, T_tank - 273.15, '--', 'LineWidth', 2.2, 'DisplayName', 'Tank water (model)');

% NEW: greenhouse air temperature (°C)
plot(t_sim, T_air  - 273.15, '-',  'LineWidth', 2.2, 'DisplayName', 'Greenhouse air (model)');

% (Optional) overlay measured inlet water from file, if you want it back:
% plot(t_txt_sec, Tin_txt, 'LineWidth', 1.8, 'DisplayName', 'Tin (measured)');

xlabel('Time (s)'); ylabel('Temperature (°C)');
title('Greenhouse model temperatures over time');
legend('Location','best'); set(gca,'FontSize',11);
