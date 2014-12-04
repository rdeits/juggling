# Juggling as a Mixed-Integer Convex Program

This package solves for a periodic trajectory of a set of N balls juggled by M hands for K time steps as a mixed-integer linear program.

## Requirements
The package requires the free Drake toolbox <https://github.com/RobotLocomotion/drake/wiki> and the Gurobi solver, which is free only for academic use <http://www.gurobi.com/>.

## Sample results

3 balls, 2 hands

```
runJuggling(3, 2, 7);
```

<a href="https://www.youtube.com/watch?v=l2Vd6WXctQc"><img src="https://raw.githubusercontent.com/rdeits/juggling/master/video/3_ball_periodic/still.png" width="75%"></a>


5 balls, 2 hands, an early solution which is highly non-optimal:

```
runJuggling(5, 2, 16);
```
<a href="https://www.youtube.com/watch?v=qz7S1uQh0Fc"><img src="https://raw.githubusercontent.com/rdeits/juggling/master/video/5_ball_screwball/still.png" width="75%"></a>

For a more optimal result, see: <https://www.youtube.com/watch?v=uOpQBBk2MJ0>.


