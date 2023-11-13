
% Simulation of CAR NK Response (No Relapse)

% Using the conditions of patient 9 (the most successful patient)
% Patient 9 was a 70 year old male, I will estimate his mass at 70kg
% He received a dose of 10m cells/kg body mass, for 700m cells total

f0=[2200.24, 0.7]; % Initial Conditions [nP0,nNK] * 10^9 Cells

rBp = 0.08; % growth rate of B-ALL cells
rNK = 0.5;  % growth rate of NKs
lNK = 0.4; %apoptosis rate of NKs
nMB = 6101.58; %carrying capacity of B-ALLs
eBP = 6; %rate of killing of B-ALLs by the NKs
KBpr = 3431.65; % Michaelis constant for effect of B-ALLs on NK growth
KBp = 7067.07; %Michaelis constant for binding of CAR to B-ALLs
KBpi = 15000; %Michaelis constant for CAR-independent binding

% Running the ode45 solver
[t,f]=ode45(@Eqs_NK_NR,0:0.1:90,f0,[], rBp, rNK, lNK, nMB, eBp, KBp, KBpr, KBpi);

LB=97.19.*f(:,1)./(1909+f(:,1)); % Leukemia tumor burden

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

% Create a new figure
figure;

yyaxis left; % Use the left y-axis
plot(t, f(:,2), 'r', 'LineWidth', 1);
ylabel('Number of NK Cells x 10^9');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

yyaxis right; % Use the right y-axis
plot(t, LB, 'b', 'LineWidth', 1);
ylabel('Tumor Burden (%)');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black
ax.YLim(1) = 0; % Set the lower limit of the y-axis to zero

title('Activated CAR NK-Cells and Tumor Burden');
xlabel('Time (days)');
legend('Activated CAR NK-Cells', 'Tumor Burden (LB)');
grid on;
hold off;

