
% Simulation of CD19- relapse.

f0=[17212.23022,0,0.3,19.89]; % Initial Conditions [nP0,nTA0,nTN0,nN0] * 10^9 Cells

Param=[0.086	
    1.7	
    0.04	
    0.073	
    19988.53	
    8.6	
    794.99	
    827.08	
    93066.8	
    0.19 
    0.1 
    0 
    7.9 
    16956.03];

rBp=Param(1);
rTA0=Param(2);
lTA0=Param(3);
lTN=Param(4); 
nMB=Param(5);
eBp=Param(6);
KBp=Param(7);
KBpr=Param(8);
KBpTN=Param(9);
ka=Param(10);    
rBn=Param(11);
km=Param(12); % mutation constant
kb=Param(13); % bystander killing
KBn=Param(14);

% all other parameters the same as in CR PR NR

[t,f]=ode45(@Eqs_NegR,[0:1:220],f0,[], rBp, rTA0, lTA0, lTN, nMB, eBp, ka, KBp, KBpr, KBpTN, rBn, km, kb, KBn);

figure;
subplot(2,2,1)
plot(t,f(:,1));
title('nB+'); % Number of B-all CD19+ cells
xlabel('Time (days)')
ylabel('Number of Cells')

subplot(2,2,2)
plot(t,f(:,2));
title('nTA'); % Number of activated CAR-T cells
xlabel('Time (days)')
ylabel('Number of Cells')

subplot(2,2,3)
plot(t,f(:,3));
title('nTN'); % Number of inactive CAR-T cells
xlabel('Time (days)')
ylabel('Number of Cells')

subplot(2,2,4)
plot(t,f(:,4));
title('nB-'); % Number of B-all CD19- cells
xlabel('Time (days)')
ylabel('Number of Cells')


LB_p=97.19.*f(:,1)./(1909+f(:,1)); % Tumor burden of B+ cells
LB_n=97.19.*f(:,4)./(1909+f(:,4)); % Tumor burden of B- cells
    
% Create a new figure for the second set of ODE solutions
figure;

yyaxis left; % Use the left y-axis
plot(t, f(:,2), 'r', 'LineWidth', 1);
ylabel('Number of Cells x 10^9');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

yyaxis right; % Use the right y-axis
plot(t, LB_p, 'b', 'LineWidth', 1);
hold on
plot(t, LB_n, '-g', 'LineWidth', 1);
ylabel('Tumor Burden (%)');

ax = gca; % Get the current axis
ax.YColor = 'k'; % Set the y-axis color to black

title('Activated CAR T-Cells and Tumor Burden');
xlabel('Time (days)');
legend('Activated CAR T-Cells', 'B+ Tumor Burden (LB_p)', 'B- Tumor Burden (LB_n)');
grid on;



