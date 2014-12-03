function obj = runJuggling(num_balls, num_hands, num_time_steps, dt)

if nargin < 4
  dt = 0.3;
end
if nargin < 3
  num_time_steps = 7;
end
if nargin < 2
  num_hands = 2;
end
if nargin < 1
  num_balls = 3;
end

checkDependency('gurobi');

obj = JugglingProblem(num_balls, num_hands, num_time_steps);
obj = obj.setup();
obj.dt = dt;
obj = obj.solveYalmip(sdpsettings('solver', 'gurobi', 'verbose', 1, 'gurobi.MIPGap', 1e-4));

h = figure(2);
obj.draw(h)

[v, xtraj] = obj.visualize();
v.playback(xtraj, struct('slider', true))

