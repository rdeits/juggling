function p = runJuggling();
checkDependency('gurobi');

s = sdpvar(2,3,4);
numel(s)

p = JugglingProblem(3, 2, 15);
p = p.setup();
p = p.solve();

h = figure(10);
p.draw(h)

