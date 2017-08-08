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

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;


// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

double car_a = 0;
double prev_car_speed = 0;

// time of prediction into the future which leads to 40 path points with t_base seconds between each
int path_vals = 50;
double t_base = 0.02;
double t_end = path_vals * t_base;

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


vector<double> JMT(vector< double> start, vector <double> end, double T)
{
  /*
  Calculate the Jerk Minimizing Trajectory that connects the initial state
  to the final state in time T.

  INPUTS

  start - the vehicles start location given as a length three array
      corresponding to initial values of [s, s_dot, s_double_dot]

  end   - the desired end state for vehicle. Like "start" this is a
      length three array.

  T     - The duration, in seconds, over which this maneuver should occur.

  OUTPUT
  an array of length 6, each value corresponding to a coefficent in the polynomial
  s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

  EXAMPLE

  > JMT( [0, 10, 0], [10, 10, 0], 1)
  [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
  */

  MatrixXd A = MatrixXd(3, 3);
  A << T*T*T, T*T*T*T, T*T*T*T*T,
       3*T*T, 4*T*T*T,5*T*T*T*T,
       6*T, 12*T*T, 20*T*T*T;

  MatrixXd B = MatrixXd(3,1);
  B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
       end[1]-(start[1]+start[2]*T),
       end[2]-start[2];

  MatrixXd Ai = A.inverse();

  MatrixXd C = Ai*B;

  vector <double> result = {start[0], start[1], .5*start[2]};
  for(int i = 0; i < C.size(); i++)
  {
    result.push_back(C.data()[i]);
  }

  return result;
}


int main() {
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
          auto car_lane = (int) floor(car_d / 3);
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

          int target_waypoint = NextWaypoint(car_x, car_y, deg2rad(car_yaw), map_waypoints_x, map_waypoints_y);

          // if waypoint is too close to current pos, choose the next one
//          if (fabs(target_x*target_x + target_y*target_y - car_x*car_x - car_y*car_y) < 250) {
//            target_waypoint = target_waypoint + 1;
//          }

          auto target_x = (double) map_waypoints_x[target_waypoint];
          auto target_y = (double) map_waypoints_y[target_waypoint];
          double target_s; double target_dx; double target_dy;
          target_dx = (double) map_waypoints_dx[target_waypoint];
          target_dy = (double) map_waypoints_dy[target_waypoint];
          // target_s = (double) map_waypoints_s[target_waypoint];
          // based on the target waypoints calculate heading based on curvature of the road
          double heading = atan2((map_waypoints_y[target_waypoint]-map_waypoints_y[target_waypoint - 1]),
                                 (map_waypoints_x[target_waypoint]-map_waypoints_x[target_waypoint - 1]));
          // calculate heading of the car from the previous path
          int path_size = previous_path_x.size(); double angle; double pos_x; double pos_y;
          if(path_size <= 1)
          {
            pos_x = car_x;
            pos_y = car_y;
            angle = deg2rad(car_yaw);
          }
          else
          {
            pos_x = (double) previous_path_x[path_size-1];
            pos_y = (double) previous_path_y[path_size-1];
            int index = 1;
            for (int i = 1; i < path_size; i++) {
              if (fabs(car_x - pos_x) > fabs(car_x - (double) previous_path_x[i])) {
                if (fabs(car_y - pos_y) > fabs(car_y - (double) previous_path_y[i])) {
                  index = i;
                  pos_x = (double) previous_path_x[index];
                  pos_y = (double) previous_path_y[index];
                }
              }
            }
            double pos_x2 = (double) previous_path_x[index - 1];
            double pos_y2 = (double) previous_path_y[index - 1];
            angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
          }


          double target_d = target_dx*target_x + target_dy*target_y;
          double car_v_s = car_speed * cos(angle);
          double car_v_d = car_speed * sin(angle);
          double car_a = (car_speed - prev_car_speed) / 0.02;
          prev_car_speed = car_speed;


          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          // search around the road in the next 5 seconds to find if a lane change is necessary

          bool change_lane = false;
          bool slow_down = false;
          double slow_down_v = 0;
          double s_to_front_car = 0;
          double obj_d;
          double obj_s;
          double obj_x;
          double obj_y;
          double obj_yaw;
          double obj_v;
          double obj_v_s;


          // if velocity is too slow, speed up in the same lane
          vector<double> XY; vector<double> SD;
          double n_s; double n_t; double x; double y; double d; double dx; double dy; double n_d;
          // change_lane = false;
          double target_v_s = 0; double target_v_d = 0;
          n_s = 0;


          for (int i = 0; i < sensor_fusion.size(); i++) {
            obj_d = sensor_fusion[i][6];
            obj_s = sensor_fusion[i][5];
            obj_x = sensor_fusion[i][1];
            obj_y = sensor_fusion[i][2];
            obj_yaw = sensor_fusion[i][4];
            // check the car is in front and in the car's lane
            if (fabs(obj_d - car_d) < 2.1 && obj_s < (car_v_s) * t_end) {
              // check if maintaining current velocity will cause a collision
              obj_yaw = sensor_fusion[i][4];
              obj_v = sensor_fusion[i][3];
              // this is an approximation, but accurate since the other cars dont seem to change lane
              if (s_to_front_car > obj_s) { s_to_front_car = obj_s; }
              slow_down = true;
              slow_down_v = obj_v * cos(deg2rad(obj_yaw));
            }
          }

          double a = 0;
          if (car_speed < 0) {
            a = -0.35;
            target_v_s = car_v_s + (a * t_end);
          } else if (car_speed < 5) {
            a = 18;
            target_v_s = car_v_s + (a * t_end);
          } else if (car_speed < 10) {
            a = 0.9;
            target_v_s = car_v_s + (a * t_end);
          } else if (car_speed < 18) {
            a = 0.35;
            target_v_s = car_v_s + (a * t_end);
          } else {
            target_v_s = 18;
            car_v_s = 18;
            car_v_d = 0;
          }
          if (slow_down) {
            target_s = s_to_front_car - 1;
            target_v_s = slow_down_v - 2;
          } else {
            target_s = car_v_s * t_end + (0.5 * a * t_end * t_end);
          }

          if (target_s < 0) {
            target_s = target_v_s * t_end;
          }

          // target_d = ((car_lane - 1) * 4) + 2.16;
          // d = target_d + 6;
          // max sure max acceleration of 2m/s**2 for s, and 0.05m/s**2 for d is complied with
          // SD = getFrenet(car_x, car_y, angle, map_waypoints_x, map_waypoints_y);

          // solve jerk minimizing trajectory for proposed goal and start states
          vector<double> jmt_s = JMT({0, car_v_s, car_a * cos(heading)}, {target_s, target_v_s, 0}, t_end);
          vector<double> jmt_d = JMT({0, car_v_d, car_a * sin(heading)}, {6.16 - car_d, 0, 0}, t_end);

          cout << 6.16 - car_d << " - " << target_s << " - " << car_x << " - " << car_y << " - " << slow_down << endl;
          cout << " ---------------------- " << endl;
          int count = 0;
          // interpolate path points
          for (int i = 0; i < path_vals; i++) {
            n_t = t_base * i;
            n_s = jmt_s[0] + jmt_s[1] * n_t + jmt_s[2] * (n_t * n_t) + jmt_s[3] * (n_t * n_t * n_t) +
                  jmt_s[4] * (n_t * n_t * n_t * n_t) + jmt_s[5] * (n_t * n_t * n_t * n_t * n_t);
            n_d = jmt_d[0] + jmt_d[1] * n_t + jmt_d[2] * (n_t * n_t) + jmt_d[3] * (n_t * n_t * n_t) +
                  jmt_d[4] * (n_t * n_t * n_t * n_t) + jmt_d[5] * (n_t * n_t * n_t * n_t * n_t);

//              dx = target_dx + (n_d * cos(a_2 - pi()/2));
//              dy = target_dy + (n_d * sin(a_2 - pi()/2));
            // the x,y,s along the segment
            // XY = getXY(car_s + n_s, car_d + n_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            x = car_x + (n_s * cos(heading)) + (n_d * cos(heading - pi()/2));
            y = car_y + (n_s * sin(heading)) + (n_d * sin(heading - pi()/2));

            cout << x << " - " << y << " - " << (n_s * cos(heading)) << " - " << n_d * cos(heading - pi()/2) << endl;
            if (count and x > next_x_vals[count]) {
              count += 1;
              next_x_vals.push_back(x);
              next_y_vals.push_back(y);
            }
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
















































































