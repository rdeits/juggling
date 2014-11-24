classdef JugglingVisualizer < Visualizer
  properties
    num_hands;
    num_balls;
    inputFrame;
    hand_colors = {'b', 'k', 'y'};
    ball_colors = {'r', 'm', 'g', 'c', [1,.8, .2]};
  end

  methods
    function obj = JugglingVisualizer(num_balls, num_hands)
      coords = {};
      for i = 1:num_balls
        coords = [coords, ...
          sprintf('ball_%d_x', i),...
          sprintf('ball_%d_y', i),...
          sprintf('ball_%d_z', i)];
      end
      for j = 1:num_hands
        coords = [coords,...
          sprintf('hand_%d_x', j),...
          sprintf('hand_%d_y', j),...
          sprintf('hand_%d_z', j),...
          sprintf('hand_%d_contact', j)];
      end
      inputFrame = CoordinateFrame('JugglingVisualizerInput', length(coords), 'x', coords);
      obj = obj@Visualizer(inputFrame);
      obj.inputFrame = inputFrame;
      obj.num_balls = num_balls;
      obj.num_hands = num_hands;
    end

    function drawWrapper(obj,t,y)
      sfigure(obj.fignum);
      draw(obj,t,y);
      if (obj.display_time)
        title(['t = ', num2str(t,'%.2f') ' sec']);
      end
      drawnow;
    end

    function draw(obj, t, x)
      persistent has_set_view
      if isempty(has_set_view)
        az = 7;
        el = 10;
        has_set_view = true;
      else
        [az, el] = view();
      end
      cla
      hold on

      p = Point(obj.inputFrame, x);

      for j = 1:obj.num_hands
        if p.(sprintf('hand_%d_contact', j))
          style = {'Color', obj.hand_colors{j}, 'MarkerSize', 15};
        else
          style = {'Color', obj.hand_colors{j}, 'MarkerSize', 10};
        end
        plot3(p.(sprintf('hand_%d_x', j)), p.(sprintf('hand_%d_y', j)), p.(sprintf('hand_%d_z', j)), 'bo', style{:});
      end

      for i = 1:obj.num_balls
        plot3(p.(sprintf('ball_%d_x', i)), p.(sprintf('ball_%d_y', i)), p.(sprintf('ball_%d_z', i)), 'ro', 'Color', obj.ball_colors{i}, 'MarkerFaceColor', obj.ball_colors{i}, 'MarkerSize', 8);
      end

      xlim([0, 1])
      ylim([-0.5, 0.5])
      zlim([-0.25, 1.5])
      % axis equal
      view(az, el);
    end
  end
end

