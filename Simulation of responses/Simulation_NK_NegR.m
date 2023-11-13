
% Simulation of CAR NK Response (No Relapse)

% Using the conditions of patient 9 (the most successful patient)
% Patient 9 was a 70 year old male, I will estimate his mass at 70kg
% He received a dose of 10m cells/kg body mass, for 700m cells total

f0=[17212.23022, 0.7, 19.89]; % Initial Conditions [nP0,nNK, nN0] * 10^9 Cells

rBp = 0.089; % growth rate of B-ALL cells
rNK = 2.00;  % growth rate of NKs
lNK = 0.08; %apoptosis rate of NKs
nMB = 19988.53; %carrying capacity of B-ALLs
eBp = 20; %rate of killing of B-ALLs by the NKs
KBpr = 1983.64; % Michaelis constant for effect of B-ALLs on NK growth
KBp = 1050.9; %Michaelis constant for binding of CAR to B-ALLs
KBpi = 10000; %Michaelis constant for CAR-independent binding
rBn = 0.1;
km = 1.5*10^-7;
kb = 17.9;
KBn = 16956.03;

% Running the ode45 solver
[t,f]=ode45(@Eqs_NK_NegR,0:0.1:90,f0,[], rBp, rNK, lNK, nMB, eBp, KBp, KBpr, KBpi, rBn, km, kb, KBn);

LB_p=97.19.*f(:,1)./(1909+f(:,1)); % Tumor burden of B+ cells
LB_n=97.19.*f(:,3)./(1909+f(:,3)); % Tumor burden of B- cells

figure;
subplot(2,2,1)
plot(t,f(:,1));
title('Number of CD19+ B-ALL Cells'); % Number of B-all leukemia cells
xlabel('Time (days)')
ylabel('Number of Cells')
hold on

subplot(2,2,2)
plot(t,f(:,2));
title('Number of CAR NK Cells'); % Number of activated CAR NK-cells
xlabel('Time (days)')
ylabel('Number of Cells')
hold off 

subplot(2,2,3)
plot(t,f(:,3));
title('nB-'); % Number of B-all CD19- cells
xlabel('Time (days)')
ylabel('Number of Cells')

% Create a new figure
figure;

yyaxis left; % Use the left y-axis
plot(t, f(:,2), 'r', 'LineWidth', 1);
ylabel('Number of NK Cells x 10^9');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

yyaxis right; % Use the right y-axis
plot(t, LB_p, 'b', 'LineWidth', 1);
ylabel('Tumor Burden (%)');
hold on
plot(t, LB_n, '-g', 'LineWidth', 1);

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black
ax.YLim(1) = 0; % Set the lower limit of the y-axis to zero

title('Activated CAR NK-Cells and Tumor Burden');
xlabel('Time (days)');
grid on;

% Add horizontal lines
yline(5, '--', 'Color', 'k', 'Label', '5%');
yline(25, '--', 'Color', 'k', 'Label', '25% (relapse threshold)');
legend('Activated CAR NK-Cells', 'B+ Tumor Burden (LB)', 'B- Tumor Burden', '5% LB', '25% LB');

hold off;

% Create a phase portrait
figure;
plot(f(:,1), f(:,2), 'LineWidth', 1);
hold on;

title('Phase Portrait with Nullclines: CD19+ B-ALL Cells vs NK Cells');
xlabel('Number of CD19+ B-ALL Cells');
ylabel('Number of NK Cells');
legend('Trajectory', 'Nullcline for nB-', 'Nullcline for nNK');
grid on;
hold off;

