classdef JugglingProblem < MixedIntegerConvexProgram
  properties
    num_hands = 2;
    num_balls = 3;
    num_frames = 15;
    dim = 3;
    dt = 0.3;
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
        obj.hand_ranges(j) = struct('center', [0.3 * j; 0; 0], 'radius', 0.1);
      end

      obj = obj.addContinuity();
      obj = obj.addHandForces();
      obj = obj.addHandRanges();
      obj = obj.addContactConstraints();
      obj = obj.addBallPlanes();
      obj = obj.addPeriodicity();

      for j = 1:obj.num_hands
        for i = 1:obj.num_balls
          obj = obj.addSymbolicConstraints([...
            1 <= sum(obj.vars.contact.symb(i,j,:)),...
            sum(obj.vars.contact.symb(i,j,:)) <= obj.num_frames-2,...
            ]);
        end
      end


      for j = 1:obj.num_hands
        % obj = obj.addSymbolicObjective(sum(sum(sum(polyderiv(obj.vars.hand_coefs.symb(:,:,:,j), 2).^2)))/1e4);
        obj = obj.addSymbolicObjective(sum(sum(sum(abs(polyderiv(obj.vars.hand_coefs.symb(:,:,:,j), 2)))))/1e4);
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

    function obj = addPeriodicity(obj)
      for j = 1:obj.num_hands
        obj = obj.addSymbolicConstraints([...
          ppval(obj.hand_traj(j), obj.breaks(1)) == ppval(obj.hand_traj(j), obj.breaks(end)),...
          ppval(fnder(obj.hand_traj(j), 1), obj.breaks(1)) == ppval(fnder(obj.hand_traj(j), 1), obj.breaks(end)),...
          ]);
      end
      for i = 1:obj.num_balls
        % next = mod(i, obj.num_balls)+1;
        next = i;
        x0 = ppval(obj.ball_traj(i), obj.breaks(1));
        x1 = ppval(obj.ball_traj(next), obj.breaks(end));
        xd0 = ppval(fnder(obj.ball_traj(i), 1), obj.breaks(1));
        xd1 = ppval(fnder(obj.ball_traj(next), 1), obj.breaks(end));
        ind = [1,2,3];
        obj = obj.addSymbolicConstraints([...
          x0(ind) == x1(ind),...
          xd0(ind) == xd1(ind),...
          ]);
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
            ];
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

    function obj = addBallPlanes(obj)
      for k = 1:obj.num_frames-1
        for i = 1:obj.num_balls
          x0 = ppval(obj.ball_traj(i), (k-1)*obj.dt);
          obj = obj.addSymbolicConstraints([...
            x0(2) == (i - obj.num_balls/2) * 0.02]);
            % x0(2) == -0.02]);
        end
      end
    end

    function trajs = extractTrajectories(obj)

      trajs = struct();

      for j = 1:obj.num_hands
        trajs.hand(j) = mkpp(obj.breaks, obj.vars.hand_coefs.value(:,:,:,j), obj.dim);
        trajs.hand_contact(j) = zoh(obj.breaks, [reshape(sum(obj.vars.contact.value(:,j,:), 1), 1, []), sum(obj.vars.contact.value(:,j,end), 1)]);
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

    function [v, xtraj] = visualize(obj)
      v = JugglingVisualizer(obj.num_balls, obj.num_hands);

      trajs = obj.extractTrajectories();
      xtraj = trajs.ball(1);
      [breaks, coefs, l, k, d] = unmkpp(trajs.ball(1));
      coefs = reshape(coefs, [d, l, k]);
      for i = 2:obj.num_balls
        [~, c2, l, k, d] = unmkpp(trajs.ball(i));
        c2 = reshape(c2, [d, l, k]);
        coefs = vertcat(coefs, c2);
      end

      for j = 1:obj.num_hands
        [~, c2, l, k, d] = unmkpp(trajs.hand(j));
        c2 = reshape(c2, [d, l, k]);
        coefs = vertcat(coefs, c2);
        [~, c2, l, k, d] = unmkpp(trajs.hand_contact(j));
        c2 = reshape(c2, [d, l, k]);
        c2 = cat(3, zeros(1, length(breaks)-1, obj.degree), c2);
        coefs = vertcat(coefs, c2);
      end

      base_coefs = coefs;
      base_breaks = breaks;
      for k = 1:15
        new_breaks = base_breaks;
        new_breaks = new_breaks + (breaks(end)-breaks(1));

        new_coefs = base_coefs;
        % for i = 1:obj.num_balls
        %   if i == 1
        %     prev = obj.num_balls;
        %   else
        %     prev = i - 1;
        %   end
        %   new_coefs(obj.dim*(i-1) + (1:obj.dim),:,:) = base_coefs(obj.dim*(prev-1) + (1:obj.dim),:,:);
        % end
        base_coefs = new_coefs;
        coefs = horzcat(coefs, new_coefs);
        breaks = horzcat(breaks, new_breaks(2:end));
      end

      xtraj = PPTrajectory(mkpp(breaks, coefs, size(coefs, 1)));
      xtraj = xtraj.setOutputFrame(v.getInputFrame());

    end

    function draw(obj, h)
      trajs = obj.extractTrajectories();
      ts = linspace(obj.breaks(1), obj.breaks(end), 100);

      hand_colors = jet(10);
      figure(2)
      clf
      subplot(211)
      hold on
      for i = 1:obj.num_balls
        x = ppval(trajs.ball(i), ts);
        plot(ts, x(3,:), 'r-');
        x = ppval(trajs.ball(i), obj.breaks);
        plot(obj.breaks, x(3,:), 'ro');
      end
      for j = 1:obj.num_hands
        x = ppval(trajs.hand(j), ts);
        plot(ts, x(3,:), 'b-', 'Color', hand_colors(j,:));
        x = ppval(trajs.hand(j), obj.breaks);
        plot(obj.breaks, x(3,:), 'b.', 'Color', hand_colors(j,:));
      end
      subplot(212)
      hold on
      for i = 1:obj.num_balls
        for j = 1:obj.num_hands
          x = ppval(trajs.hand_ball_force(i,j), ts);
          plot(ts, x(3,:), 'b-');
        end
      end
    end
  end
end
