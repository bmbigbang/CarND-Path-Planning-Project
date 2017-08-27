# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program
   
Using a finite state machine to control the car's intentions for the next 1
second, we can track the car's driving behaviour and adjust the path points given to the car
to obey the rules of the highway and avoid collisions. At each point of time the program
explores the car's nearby objects given by the sensor fusion array. Regardless of the state of the car
the reference speed is reduced if an object is observed in front of the car (main.cpp lines 364-367).
The objects in the other lanes are also tracked and based on their distance from
the car, it is decided whether that lane is a safe lane to change into (main.cpp lines 382-385).

Once the state of the system is known, a hierarchy of values is used to manage and
maintain the speed of the car (main.cpp lines 419-424). The speed of the car in turn is 
used to calculate the 
trajectory, which is fed back to the simulator on a 0.02 time difference basis. Therefore
this speed is the defining factor in the car's acceleration. The spline library is used
to find closest points that would give us a smooth line following the waypoints given
by the map. Depending on the state of the car found earlier, we adjust the reference speed
of the car to control its speed along the s direction. For the d direction, simply the target
lane value is given to the spline library which can cause a big jump in the car's motion if a big 
enough difference is given. As a result it is attempted to increment this change in smaller sections.
Once a lane change state has been initialised, the d value of the spline is set to the new lane minus some small
increment (main.cpp line 431). Once the car goes to the edge of the lane, a lane change is found and the car 
returns to the normal state. In this state the car is constantly correcting its d value every second 
(main.cpp lines 444-451). This is necessary because the waypoints are not fine enough to resolve the entire 
highway and at some corners will force us out of the lane slightly.

The spline tool is then used, alongside an estimate of the distance between each path point in the future 
(main.cpp lines 491-494), to calculate a close estimate of future path points y values. The previous 
path points that are unused will also be sent back for the future path points to ensure a continuous path
is given to the car at each state. This also means predicting differences in less than 1 second in the future will 
not always be possible.
   
### Simulator. You can download the Term3 Simulator BETA which contains the Path Planning Project from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.



## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

