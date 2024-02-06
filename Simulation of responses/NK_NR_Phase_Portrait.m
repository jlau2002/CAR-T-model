function NK_NR_Phase_Portrait

% Define parameters (replace with your actual parameters)

rBp = 0.08; % growth rate of B-ALL cells
rNK = 0.5;  % growth rate of NKs
lNK = 0.4; %apoptosis rate of NKs
nMB = 6101.58; %carrying capacity of B-ALLs
eBP = 6; %rate of killing of B-ALLs by the NKs
KBpr = 3431.65; % Michaelis constant for effect of B-ALLs on NK growth
KBp = 7067.07; %Michaelis constant for binding of CAR to B-ALLs
KBpi = 15000; %Michaelis constant for CAR-independent binding


% Fix initial number of NK cells
initial_nNK = 0.7;

% Generate 20 different initial conditions for nB
initial_nB_range = rand(20, 1) * 10000 + 5000; % You can adjust the range based on your system

% Create a figure for the phase portrait
figure;
hold on;

% Iterate through different initial conditions for nB
for i = 1:size(initial_nB_range, 1)
    % Run the ode solver
    [t, f] = ode45(@Eqs_NK_CR, 0:0.1:1000, [initial_nB_range(i), initial_nNK], [], rBp, rNK, lNK, nMB, eBP, KBp, KBpr, KBpi);

    % Plot the trajectory
    plot(f(:, 2), f(:, 1), 'LineWidth', 1);

    hold on;
    
end

% Set y-axis minimum value to zero
ylim([0, inf]);

title('Phase Portrait: CD19+ B-ALL Cells vs NK Cells');
xlabel('Number of NK Cells * 10^9');
ylabel('Number of B-ALL Cells * 10^9');
grid on;
hold off;
end