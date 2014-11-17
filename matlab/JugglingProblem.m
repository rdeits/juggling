classdef JugglingProblem < MixedIntegerConvexProgram
  properties
    num_hands = 2;
    num_balls = 3;
    num_frames = 15;
    dim = 3;
    dt = 0.1;
    hand_traj = struct('form', {}, 'breaks', {}, 'coefs', {}, 'pieces', {}, 'order', {}, 'dim', {});
    ball_traj = struct('form', {}, 'breaks', {}, 'coefs', {}, 'pieces', {}, 'order', {}, 'dim', {});
    % hand_force_traj;
    hand_ranges;
    breaks;
    degree = 3;
    g = 9.81;
    ball_mass = 0.1;
  end

  methods
    function obj = JugglingProblem(num_balls, num_hands, num_frames)
      obj = obj@MixedIntegerConvexProgram(true);
      if ~isempty(num_balls)
        obj.num_balls = num_balls;
      end
      if ~isempty(num_hands)
        obj.num_hands = num_hands;
      end
      if ~isempty(num_frames)
        obj.num_frames = num_frames;
      end
    end

    function obj = setup(obj)
      obj.breaks = 0:obj.dt:(obj.dt*(obj.num_frames-1));
      obj = obj.addVariable('hand_coefs', 'C', [obj.dim, obj.num_frames-1, obj.degree+1, obj.num_hands], -100, 100);
      obj = obj.addVariable('ball_coefs', 'C', [obj.dim, obj.num_frames-1, obj.degree+1, obj.num_balls], -100, 100);

      for j = 1:obj.num_hands
        obj.hand_traj(j) = mkpp(obj.breaks,obj.vars.hand_coefs.symb(:,:,:,j), obj.dim);
      end 
      for j = 1:obj.num_balls
        obj.ball_traj(j) = mkpp(obj.breaks,obj.vars.ball_coefs.symb(:,:,:,j), obj.dim);
      end

      obj.hand_ranges = struct('center', {}, 'radius', {});
      for j = 1:obj.num_hands
        obj.hand_ranges(j) = struct('center', [0; 0.3 * j; 0], 'radius', 0.1);
      end


      obj = obj.addContinuity();
      obj = obj.addHandForces();
      obj = obj.addHandRanges();
      obj = obj.addContactConstraints();
      obj = obj.addInitialState();
      obj = obj.addBallRange();

      for i = 1:obj.num_balls
        x0 = ppval(obj.ball_traj(i), 0);
        obj.symbolic_constraints = [obj.symbolic_constraints,...
          x0(3) <= 1,...
          % ppval(obj.ball_traj, 0) == ppval(obj.hand_traj, 0),...
          ppval(fnder(obj.ball_traj(i), 1), 0) == 0,...
          ];
      end

      for j = 1:obj.num_hands
        obj = obj.addSymbolicObjective(sum(sum(squeeze(polyderiv(obj.vars.hand_coefs.symb, 3)).^2, 1)));
      end

    end
    
    function obj = addContinuity(obj)
      for k = 1:obj.num_frames-2
        for dorder = 0:1
          for j = 1:obj.num_hands
            c1 = reshape(polyderiv(obj.vars.hand_coefs.symb(:,k,:,j), dorder), obj.dim, []);
            c2 = reshape(polyderiv(obj.vars.hand_coefs.symb(:,k+1,:,j), dorder), obj.dim, []);
            obj.symbolic_constraints = [obj.symbolic_constraints,...
              c1 * (obj.dt.^(obj.degree-dorder:-1:0))' == c2 * (0.^(obj.degree-dorder:-1:0))',...
              ];
          end
          for i = 1:obj.num_balls
            c1 = reshape(polyderiv(obj.vars.ball_coefs.symb(:,k,:,i), dorder), obj.dim, []);
            c2 = reshape(polyderiv(obj.vars.ball_coefs.symb(:,k+1,:,i), dorder), obj.dim, []);
            obj.symbolic_constraints = [obj.symbolic_constraints,...
              c1 * (obj.dt.^(obj.degree-dorder:-1:0))' == c2 * (0.^(obj.degree-dorder:-1:0))',...
              ];
          end
        end
      end
    end

    function obj = addHandRanges(obj)
      for j = 1:obj.num_hands
        for k = 1:obj.num_frames-1
          obj.symbolic_constraints = [obj.symbolic_constraints,...
            abs(ppval(obj.hand_traj(j), (k-1)*obj.dt) - obj.hand_ranges(j).center) <= obj.hand_ranges(j).radius,...
            abs(ppval(obj.hand_traj(j), (k)*obj.dt) - obj.hand_ranges(j).center) <= obj.hand_ranges(j).radius,...
            ];
        end
      end
    end

    function obj = addHandForces(obj)
      grav_coefs = zeros(obj.dim, obj.degree-1);
      grav_coefs(end, end) = -obj.g;

      obj = obj.addVariable('hand_ball_force_coefs', 'C', [obj.dim, obj.num_frames-1, obj.degree-1, obj.num_balls, obj.num_hands], -100, 100);

      for k = 1:obj.num_frames-1
        for i = 1:obj.num_balls
          total_force_coefs = reshape(sum(obj.vars.hand_ball_force_coefs.symb(:,k,:,i,:), 5), obj.dim, []);
          ball_accel_coefs = reshape(polyderiv(obj.vars.ball_coefs.symb(:,k,:,i), 2), obj.dim, []);
          total_force_coefs = total_force_coefs + obj.ball_mass * grav_coefs;
          obj.symbolic_constraints = [obj.symbolic_constraints,...
            obj.ball_mass * ball_accel_coefs == total_force_coefs,...
            total_force_coefs(1:2,:) == 0,...
            ];
        end
      end
    end

    function obj = addInitialState(obj)
      for j = 1:obj.num_hands
        for k = 1:obj.num_frames-1
          obj.symbolic_constraints = [obj.symbolic_constraints,...
            ppval(obj.hand_traj(j), 0) == obj.hand_ranges(j).center];
        end
      end
    end

    function obj = addContactConstraints(obj)
      obj = obj.addVariable('contact', 'B', [obj.num_balls, obj.num_hands, obj.num_frames-1], 0, 1);
      contact = obj.vars.contact.symb;
      for k = 1:obj.num_frames-1
        obj.symbolic_constraints = [obj.symbolic_constraints,...
          sum(contact(:,:,k), 1) <= 1,...
          ];
        for j = 1:obj.num_hands
          for i = 1:obj.num_balls
            obj.symbolic_constraints = [obj.symbolic_constraints,...
              implies(~contact(i,j,k), obj.vars.hand_ball_force_coefs.symb(:,k,:,i,j) == 0),...
              implies(contact(i,j,k), obj.vars.hand_coefs.symb(:,k,:,j) == obj.vars.ball_coefs.symb(:,k,:,i)),...
              ];
          end
        end
      end
    end

    function obj = addBallRange(obj)
      for k = 1:obj.num_frames-1
        for i = 1:obj.num_balls
          x0 = ppval(obj.ball_traj(i), (k-1)*obj.dt);
          x1 = ppval(obj.ball_traj(i), (k)*obj.dt);
          obj.symbolic_constraints = [obj.symbolic_constraints,...
            [-1; -1; 0] <= x0,...
            x0 <= [2; 1; inf],...
            [-1; -1; 0] <= x1,...
            x1 <= [2; 1; inf],...
            ];
        end
      end
    end

    function trajs = extractTrajectories(obj)

      trajs = struct();

      for j = 1:obj.num_hands
        trajs.hand(j) = mkpp(obj.breaks, obj.vars.hand_coefs.value(:,:,:,j), obj.dim);
      end
      for j = 1:obj.num_balls
        trajs.ball(j) = mkpp(obj.breaks, obj.vars.ball_coefs.value(:,:,:,j), obj.dim);
      end
      for j = 1:obj.num_hands
        for i = 1:obj.num_balls
          trajs.hand_ball_force(i,j) = mkpp(obj.breaks, obj.vars.hand_ball_force_coefs.value(:,:,:,i,j), obj.dim);
        end
      end
    end

    function draw(obj, h)
      figure(h);
      clf;
      hold on
      trajs = obj.extractTrajectories();
      ts = linspace(obj.breaks(1), obj.breaks(end), 100);
      for j = 1:obj.num_hands
        x = ppval(trajs.hand(j), ts);
        plot3(x(1,:), x(2,:), x(3,:), 'b.-');
        % plot(ts, ppval(trajs.hand(j), ts), 'b.-');
      end
      for j = 1:obj.num_balls
        x = ppval(trajs.ball(j), ts);
        plot3(x(1,:), x(2,:), x(3,:), 'ro-');
        % plot(ts, ppval(trajs.ball(j), ts), 'ro-');
      end

      figure(2)
      clf
      subplot(211)
      hold on
      for i = 1:obj.num_balls
        x = ppval(trajs.ball(i), ts);
        plot(ts, x(3,:), 'b.-');
      end
      x = ppval(trajs.hand(1), ts);
      plot(ts, x(3,:), 'r.-');
      subplot(212)
      hold on
      for i = 1:obj.num_balls
        x = ppval(trajs.hand_ball_force(i,1), ts);
        plot(ts, x(3,:), 'b.-');
      end
    end
  end
end
