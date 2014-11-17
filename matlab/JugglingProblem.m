classdef JugglingProblem < MixedIntegerConvexProgram
  properties
    num_hands = 2;
    num_balls = 3;
    num_frames = 15;
    dim = 3;
    dt = 0.1;
    % hand_traj = struct('form', {}, 'breaks', {}, 'coefs', {}, 'pieces', {}, 'order', {}, 'dim', {});
    % ball_traj = struct('form', {}, 'breaks', {}, 'coefs', {}, 'pieces', {}, 'order', {}, 'dim', {});
    % hand_force_traj;
    hand_ranges;
    breaks;
    degree = 3;
    g = 9.81;
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

      % for j = 1:obj.num_hands
      %   obj = obj.addVariable(sprintf('hand_%d_coefs', j), 'C', [obj.dim, obj.num_frames-1, obj.degree+1], 0, 1);
      %   obj.hand_traj(j) = mkpp(obj.breaks,obj.vars.(sprintf('hand_%d_coefs', j)).symb, obj.dim);
      % end 
      % for j = 1:obj.num_balls
      %   obj = obj.addVariable(sprintf('ball_%d_coefs', j), 'C', [obj.dim, obj.num_frames-1, obj.degree+1], 0, 1);
      %   obj.ball_traj(j) = mkpp(obj.breaks,obj.vars.(sprintf('ball_%d_coefs', j)).symb, obj.dim);
      % end

      obj = obj.addHandForces();
      obj = obj.addHandRanges();
      obj = obj.addContactConstraints();
      obj = obj.addInitialState();
    end

    function obj = addHandRanges(obj)
      obj.hand_ranges = struct('center', {}, 'radius', {});
      for j = 1:obj.num_hands
        obj.hand_ranges(j) = struct('center', [0; 0.3 * j; 0], 'radius', 0.1);
        for k = 1:obj.num_frames-1
          obj.symbolic_constraints = [obj.symbolic_constraints,...
            % abs(ppval(obj.hand_traj(j), (k-0.5)*obj.dt) - obj.hand_ranges(j).center) <= obj.hand_ranges(j).radius];
            % ppval(obj.hand_traj(j), k*obj.dt) 
            abs(squeeze(obj.vars.hand_coefs.symb(:,k,:,j)) * [0; 0; 1] - obj.hand_ranges(j).center) <= obj.hand_ranges(j).radius];
        end
      end
    end

    function obj = addHandForces(obj)
      obj.hand_force_traj = {};
      grav_coefs = zeros(obj.dim, obj.num_frames-1, obj.degree-1);
      grav_coefs(end, :, end) = -obj.g;

      for j = 1:obj.num_hands
        obj.hand_force_traj{j} = {};
        for i = 1:obj.num_balls
          name = sprintf('hand_%d_ball_%d_force_coefs', j, i);
          obj = obj.addVariable(name, 'C', [obj.dim, obj.num_frames-1, obj.degree-1], -100, 100);
          obj.hand_force_traj{j}{i} = mkpp(obj.breaks, obj.vars.(name).symb, obj.dim);
        end
      end
    end

    function obj = addInitialState(obj)
      for j = 1:obj.num_hands
        obj.symbolic_constraints = [obj.symbolic_constraints,...
          ppval(obj.hand_traj(j), 0) == obj.hand_ranges(j).center];
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
            coefs = obj.vars.(sprintf('hand_%d_ball_%d_force_coefs', j, i)).symb;
            obj.symbolic_constraints = [obj.symbolic_constraints,...
              implies(~contact(i, j, k), coefs == 0);
              implies(contact(i, j, k), obj.vars.(sprintf('hand_%d_coefs', j)).symb == obj.vars.(sprintf('ball_%d_coefs', i)).symb),...
              ];
          end
        end
      end
    end


    function trajs = extractTrajectories(obj)

      trajs = struct();

      for j = 1:obj.num_hands
        trajs.hand(j) = mkpp(obj.breaks, obj.vars.(sprintf('hand_%d_coefs', j)).value, 3);
      end
      for j = 1:obj.num_balls
        trajs.ball(j) = mkpp(obj.breaks, obj.vars.(sprintf('ball_%d_coefs', j)).value, 3);
      end
    end

    function draw(obj, h)
      figure(h);
      clf;
      hold on
      trajs = obj.extractTrajectories();
      ts = linspace(0, obj.dt * obj.num_frames-1, 100);
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
    end
  end
end
