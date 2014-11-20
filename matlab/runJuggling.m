function obj = runJuggling();
checkDependency('gurobi');

s = sdpvar(2,3,4);
numel(s)

obj = JugglingProblem(3, 2, 9);
obj = obj.setup();
obj = obj.solve();

h = figure(2);
obj.draw(h)

[v, xtraj] = obj.visualize();
v.playback(xtraj, struct('slider', true))

