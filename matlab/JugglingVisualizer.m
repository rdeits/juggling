classdef JugglingVisualizer < Visualizer
  properties
    num_hands;
    num_balls;
    inputFrame;
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
      view(7, 10);
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
      % [az, el] = view();
      az = 7;
      el = 10;
      cla
      hold on

      p = Point(obj.inputFrame, x);

      for i = 1:obj.num_balls
        plot3(p.(sprintf('ball_%d_x', i)), p.(sprintf('ball_%d_y', i)), p.(sprintf('ball_%d_z', i)), 'ro');
      end

      for j = 1:obj.num_hands
        if p.(sprintf('hand_%d_contact', j))
          style = {'MarkerFaceColor', 'b'};
        else
          style = {};
        end
        plot3(p.(sprintf('hand_%d_x', j)), p.(sprintf('hand_%d_y', j)), p.(sprintf('hand_%d_z', j)), 'bo', style{:});
      end
      xlim([0, 1])
      ylim([-0.5, 0.5])
      zlim([-0.25, 1.5])
      % axis equal
      view(az, el);
    end
  end
end

