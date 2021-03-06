#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "json.hpp"
#include <time.h>
#include "spline.h"

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;


// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// define variables for tracking car's acceleration, time and the previous path size
double car_a = 0;
double prev_car_speed = 0;
double prev_prev_size = 0;
time_t now;
time_t t;

double seconds;
int target_lane = 1;
int new_lane = 1;

double ref_vel = 0;

// time of prediction into the future which leads to 50 path points with t_base seconds between each
int path_vals = 50;
double t_base = 0.02;
double t_end = path_vals * t_base;
bool change_lane = false;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


int main() {
  time(&t);
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];
          int car_lane = car_d / 4;

          /// target velocity should be 25 m/s for s and small value for d
          double target_s_v = 25;


          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];
          // cout << car_x << "," << car_y << "," << car_s << "," << car_d << "," << car_yaw << "," << car_speed << "," << endl;
          //for (int i = 0; i < sensor_fusion.size(); i++) {
            // cout << sensor_fusion[i] << "," << sensor_fusion.size() << endl;
            // objects orientation: [ID, x, y, speed, yaw, d, s]
          //}

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // calculate heading of the car from the previous path
          int prev_size = previous_path_x.size();
          vector<double> ptsx;
          vector<double> ptsy;

          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          // reuse previous path points where possible to deduce the yaw angle
          // and pass on points to be used in the future
          if (prev_size < 2) {
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);

            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          }
          else {
            ref_x = previous_path_x[prev_size - 1];
            ref_y = previous_path_y[prev_size - 1];

            double ref_x_prev = previous_path_x[prev_size - 2];
            double ref_y_prev = previous_path_y[prev_size - 2];
            ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);

            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);
          }

          // deduce the car acceleration and component of velocity in the s/d directions
          double car_v_s = car_speed * cos(ref_yaw);
          double car_v_d = car_speed * sin(ref_yaw);
          if (prev_prev_size - prev_size !=  0) {
            // delta v over delta t
            car_a = (car_speed - prev_car_speed) / ((prev_prev_size - prev_size) * 0.02);
          }
          else {
            car_a = 0;
          }
          prev_car_speed = car_speed;

          // define variables to use for sensor fusion and the states of the car
          bool slow_down = false;
          double obj_d;
          double obj_s;
          double obj_x;
          double obj_y;
          double obj_vx;
          double obj_vy;
          double obj_speed;
          double obj_v_s;
          int obj_lane;
          int first_obj_id = 0;
          double a = 2.5;

          // keep track of available lanes for lange changes. disallow double lane shifts if we are on an edge lane
          vector<int> available_lanes = {0, 1, 2};
          if ((car_lane == 0) || (car_lane == 2)) {
            available_lanes.erase(available_lanes.begin() + 2);
            available_lanes.erase(available_lanes.begin());
          }
          else {
            available_lanes.erase(available_lanes.begin() + car_lane);
          }

          // search around the road in the next 1 second to find if a slow down or lane change is necessary
          for (int i = 0; i < sensor_fusion.size(); i++) {
            obj_d = sensor_fusion[i][6];
            // if object is on our side of the road
            if (obj_d > 0) {
              obj_s = sensor_fusion[i][5];
              obj_lane = obj_d / 4;
              obj_vx = sensor_fusion[i][3];
              obj_vy = sensor_fusion[i][4];
              obj_speed = sqrt(obj_vx * obj_vx + obj_vy * obj_vy);
              obj_s += ((double) prev_size * .02 * obj_speed);

              // check the car is in front and in the car's lane
              if (obj_lane == car_lane) {

                if ((obj_s > car_s) && (5 < (obj_s - car_s)) && ((obj_s - car_s) < 45)) {
                  slow_down = true;
                }
              }
              else {
                // watch for swerving cars into our lane
                if ((obj_speed > ref_vel + 10) && ((car_d + 2 < obj_d) && (obj_d > car_d - 2)) &&
                    (obj_s > car_s) && (20 < (obj_s - car_s)) && ((obj_s - car_s) < 40)) {
                  slow_down = true;
                  for (int k = 0; k < available_lanes.size(); k++) {
                    if (available_lanes[k] == obj_lane) {
                      available_lanes.erase(available_lanes.begin() + k);
                    }
                  }
                }
                // look for clear lanes to change lane into
                if ((0 < (obj_s - car_s)) && ((obj_s - car_s) < 50)) {
                  for (int k = 0; k < available_lanes.size(); k++) {
                    if (available_lanes[k] == obj_lane) {
                      available_lanes.erase(available_lanes.begin() + k);
                    }
                  }
                }
              }
            }
          }

          // if there are available lanes and we are slowing down, and the car is in the middle of its own lane
          // and we are moving slowly enough, enter the change lane state
          if (!change_lane && slow_down && !available_lanes.empty() &&
              !((1 < (int(car_d) % 4)) && ((int(car_d) % 4) < 3)) && car_speed < 36) {
            time(&now);
            double last_time = difftime(now,t);
            // do not keep changing lanes in smaller windows than 5 seconds
            if (last_time > 5) {
              change_lane = true;
              slow_down = true;
              new_lane = available_lanes[0];
              target_lane = new_lane;
              time(&t);
            }
          }

          // if we are in the change lane state, correct the new lane
          if (change_lane && (target_lane != car_lane)) {
            new_lane = target_lane;
          }
          // if we are in the change lane state but have moved to the target lane, return to normal state
          else if (change_lane && (target_lane == car_lane)) {
            change_lane = false;
          }

          // hierarchy of slowing down per speed ranges and depending on the current state
          if (slow_down && change_lane && ref_vel > 30) { ref_vel -= .284; }
          else if (slow_down && ref_vel > 45) { ref_vel -= .224; }
          else if (slow_down && ref_vel > 36) { ref_vel -= .164; }
          else if (slow_down && ref_vel > 30) { ref_vel -= .134; }
          else if (ref_vel < 30) { ref_vel += 0.33; }
          else if (ref_vel < 45) { ref_vel += .224; }

          // add waypoints to our estimated trajectory spline depending on whether we are changing lanes or not
          double future_d_1 = 2 + (4 * car_lane);
          double future_d_2 = 2 + (4 * car_lane);
          double future_d_3 = 2 + (4 * car_lane);
          if (new_lane != car_lane) {
            if ((2 > (car_d - (new_lane * 4 + 2.08)) && (car_d - (new_lane * 4 + 2.08)) > -2)) {
              future_d_1 = 2 + (4 * new_lane);
              future_d_2 = 2 + (4 * new_lane);
              future_d_3 = 2 + (4 * new_lane);
            }
            else {
              // smaller increment to the current lane to allow smoother lane transition
              future_d_1 = 2 + (4 * new_lane) - (1.18 * (new_lane - car_lane));
              future_d_2 = 2 + (4 * new_lane) - (0.5 * (new_lane - car_lane));
              future_d_3 = 2 + (4 * new_lane);
            }
          }
          // if not in the changing lane position, correct the target_d to be further in the middle of the current lane
          else if (1.9 < (car_d - (car_lane * 4 + 2))) {
            future_d_1 -= (car_d - (car_lane * 4 + 2)) / 3;
            future_d_2 -= (car_d - (car_lane * 4 + 2)) / 3;
          }
          else if (-1.9 < (car_d - (car_lane * 4 + 2))) {
            future_d_1 += (car_d - (car_lane * 4 + 2)) / 3;
            future_d_2 += (car_d - (car_lane * 4+ 2)) / 3;
          }

          vector<double> next_wp0 = getXY(car_s + 30, future_d_1, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s + 60, future_d_2, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s + 90, future_d_3, map_waypoints_s, map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);

          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          for (int i = 0; i < ptsx.size(); i++) {
            double shift_x = ptsx[i] - ref_x;
            double shift_y = ptsy[i] - ref_y;
            // rotate each point in accordance to the car's orientation such that x points along the direction of s
            ptsx[i] = (shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw));
            ptsy[i] = (shift_x * sin(-ref_yaw) + shift_y * cos(-ref_yaw));
          }

          // initialise the spline fitter and add points
          tk::spline spl;
          spl.set_points(ptsx, ptsy);

          // reuse the old points for the future path
          for (int i = 0; i < previous_path_x.size(); i++) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          double target_x = 30.0;
          double target_y = spl(target_x);
          double target_dist = sqrt((target_x * target_x) + (target_y * target_y));

          double x_add_on = 0;
          double y_point;
          prev_prev_size = prev_size;
          // find and calculate the future path points in accordance to the current spline
          for (int i = 0; i <= (50 - prev_size); i++) {
            double N = target_dist / (.02 * ref_vel / 2.24);
            double x_point = x_add_on + (target_x / N);
            double y_point = spl(x_point);

            x_add_on = x_point;

            double x_ref = x_point;
            double y_ref = y_point;

            x_point = (x_ref * cos(ref_yaw)) - (y_ref * sin(ref_yaw));
            y_point = (x_ref * sin(ref_yaw)) + (y_ref * cos(ref_yaw));

            x_point += ref_x;
            y_point += ref_y;

            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);
          }


          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          // this_thread::sleep_for(chrono::milliseconds(160));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































