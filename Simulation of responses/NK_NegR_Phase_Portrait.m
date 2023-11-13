function NK_NegR_Phase_Portrait

% Define parameters
rBp = 0.089; % growth rate of B-ALL cells
rNK = 2.00;  % growth rate of NKs
lNK = 0.08; % apoptosis rate of NKs
nMB = 19988.53; % carrying capacity of B-ALLs
eBp = 20; % rate of killing of B-ALLs by the NKs
KBpr = 1983.64; % Michaelis constant for effect of B-ALLs on NK growth
KBp = 1050.9; % Michaelis constant for binding of CAR to B-ALLs
KBpi = 10000; % Michaelis constant for CAR-independent binding
rBn = 0.1;
km = 1.5 * 10^-7;
kb = 17.9;
KBn = 16956.03;

% Create a figure for the 3D phase portrait
figure;
hold on;

% Generate a grid of initial conditions
nB_range = linspace(5000, 15000, 2); % You can adjust the range based on your system
nN_range = linspace(100, 1000, 2); % You can adjust the range based on your system

% Solve ODEs and plot trajectory in 3D space
for i = 1:length(nB_range)
    for j = 1:length(nN_range)
        % Run the ode solver
        [t, f] = ode45(@Eqs_NK_NegR, 0:0.1:50, [nB_range(i), 0, nN_range(j)], [], rBp, rNK, lNK, nMB, eBp, KBp, KBpr, KBpi, rBn, km, kb, KBn);
        x=f(:,1);
        y=f(:,3);
        tbl = table(x,y,t);

        % Plot the trajectory in 3D
        plot3(tbl, "x", "y", "t", 'LineWidth', 1);
    end
end

title('3D Phase Portrait: CD19+ B-ALL Cells, CD19- B-ALL Cells, and NK Cells');
xlabel('CD19+ B-ALL Cells');
ylabel('CD19- B-ALL Cells');
zlabel('NK Cells');
grid on;
hold off;

end
