function obj = testDynamics()

num_frames = 4;
dt = 1;
obj = MixedIntegerConvexProgram(true);
obj = obj.addVariable('force_coefs', 'C', [3, num_frames-1, 2], -100, 100);
obj = obj.addVariable('ball_coefs', 'C', [3, num_frames-1, 4], -100, 100);

grav_coefs = zeros(3, 1, 3-1);
grav_coefs(end,:, end) = -9.81;

for k = 1:num_frames-1
  total_force_coefs = obj.vars.force_coefs.symb(:,k,:) + grav_coefs;
  accel_coefs = polyderiv(obj.vars.ball_coefs.symb(:,k,:), 2);
  obj = obj.addSymbolicConstraints([accel_coefs == total_force_coefs]);
end

for k = 1:num_frames-2
  for dorder = 0:1
    c1 = reshape(polyderiv(obj.vars.ball_coefs.symb(:,k,:), dorder), 3, []);
    c2 = reshape(polyderiv(obj.vars.ball_coefs.symb(:,k+1,:), dorder), 3, []);
    obj.symbolic_constraints = [obj.symbolic_constraints,...
      c1 * (dt.^(3-dorder:-1:0))' == c2 * (0.^(3-dorder:-1:0))',...
      ];
  end
end

obj = obj.addSymbolicObjective(sum(sum(polyderiv(obj.vars.ball_coefs.symb, 3).^2, 1)));


obj = obj.solve();
breaks = 0:num_frames-1;
force_traj = mkpp(breaks, obj.vars.force_coefs.value, 3);
ball_traj = mkpp(breaks, obj.vars.ball_coefs.value, 3);

obj.vars.force_coefs.value
obj.vars.ball_coefs.value

ts = linspace(breaks(1), breaks(end));
x = ppval(ball_traj, ts);
xdd = ppval(fnder(ball_traj, 2), ts);
forces = bsxfun(@plus, ppval(force_traj, ts), [0;0;-9.81]);
valuecheck(xdd, forces);

clf
subplot(211)
plot(ts, x(3,:));
subplot(212)
plot(ts, forces(3,:))