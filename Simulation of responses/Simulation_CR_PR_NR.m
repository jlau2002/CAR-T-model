
% Simulation of continuous remission, CD19+ relapse and non-response

f0=[2200.24,0,16.46]; % Initial Conditions [nP0,nTA,nT0] * 10^9 Cells

rBp=0.069; % Growth rate of B+ cells
rTA0=1.62; % Growth rate of activated T cells
lTA0=0.12; % Apoptosis rate of activated T cells
lTN=0.00003; % Apoptosis rate of inactivated T cells
nMB=2939.1; % Called "n_C" in the paper, don't know why they changed the name
eBp=22.72; % Killing rate of CD19+ B-ALLs by the CAR T-cells
ka=0.65; % Rate of conversion of inactive CAR-Ts into activated CAR-Ts
KBp=5891.45; % Michaelis constant for the killing of B-all cells and saturation of B-ALLs
KBpr=637.64; % Michaelis constant for effect of B-ALLs on growth of the CAR-Ts (called "K_r" in paper)
KBpTN=1808.02; % Michaelis constant for effect of B-ALLs on activation of the inactive CAR-Ts

[t,f]=ode45(@Eqs_CR_PR_NR,0:0.1:90,f0,[], rBp, rTA0, lTA0, lTN, nMB, eBp, ka, KBp, KBpr, KBpTN);

%figure;
subplot(2,2,1)
plot(t,f(:,1));
title('Number of CD19+ B-ALL Cells'); % Number of B-all leukemia cells
xlabel('Time (days)')
ylabel('Number of Cells')
hold on

subplot(2,2,2)
plot(t,f(:,2));
title('Number of Activated CAR T-Cells'); % Number of activated CAR T-cells
xlabel('Time (days)')
ylabel('Number of Cells')
hold on

subplot(2,2,3)
plot(t,f(:,3));
title('Number of Inactive CAR T-Cells'); % Number of inactive CAR T-cells
xlabel('Time (days)')
ylabel('Number of Cells')

hold on

LB=97.19.*f(:,1)./(1909+f(:,1)); % Leukemia tumor burden

hold off

% Create a new figure
figure;

yyaxis left; % Use the left y-axis
plot(t, f(:,2), 'r', 'LineWidth', 1);
ylabel('Number of Cells x 10^9');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

yyaxis right; % Use the right y-axis
plot(t, LB, 'b', 'LineWidth', 1);
ylabel('Tumor Burden (%)');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

title('Activated CAR T-Cells and Tumor Burden');
xlabel('Time (days)');
legend('Activated CAR T-Cells', 'Tumor Burden (LB)');
grid on;
hold off;

% Create a new figure
figure;

yyaxis left; % Use the left y-axis
semilogy(t, f(:,2), 'r', 'LineWidth', 1);
ylabel('Number of Cells x 10^9');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

yyaxis right; % Use the right y-axis
semilogy(t, LB, 'b', 'LineWidth', 1);
ylabel('Tumor Burden (%)');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

title('Activated CAR T-Cells and Tumor Burden');
xlabel('Time (days)');
legend('Activated CAR T-Cells', 'Tumor Burden (LB)');
grid on;




