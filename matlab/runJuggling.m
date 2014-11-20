function obj = runJuggling(num_balls, num_hands, num_time_steps);
checkDependency('gurobi');

s = sdpvar(2, 1,6);
numel(s)

obj = JugglingProblem(num_balls, num_hands, num_time_steps);
obj = obj.setup();
obj = obj.solve();

h = figure(2);
obj.draw(h)

[v, xtraj] = obj.visualize();
v.playback(xtraj, struct('slider', true))

