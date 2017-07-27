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
          double car_v_s = sin(car_yaw) * car_speed;
          double car_v_d = cos(car_yaw) * car_speed;
          /// target velocity should be 25 m/s for s and small value for d
          double target_s_v = 25;
          double target_d_v = 0.5;
          // time of prediction into the future which leads to 40 path points with 0.02s between each
          int path_vals = 10;
          double t_end = path_vals * 0.02;

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

          vector<double> next_x_vals(path_vals, 0);
          vector<double> next_y_vals(path_vals, 0);

          int target_waypoint = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);

          auto target_x = (double) map_waypoints_x[target_waypoint];
          auto target_y = (double) map_waypoints_y[target_waypoint];
          double target_s; double target_dx; double target_dy;

          // if waypoint is too close to current pos, choose the next one
          if (target_x*target_x + target_y*target_y - car_x*car_x - car_y*car_y < 300) {
            target_x = (double) map_waypoints_x[target_waypoint + 1];
            target_y = (double) map_waypoints_y[target_waypoint + 1];
            target_dx = (double) map_waypoints_dx[target_waypoint + 1];
            target_dy = (double) map_waypoints_dy[target_waypoint + 1];
            target_s = (double) map_waypoints_s[target_waypoint + 1];
          }
          else {
            target_dx = (double) map_waypoints_dx[target_waypoint];
            target_dy = (double) map_waypoints_dy[target_waypoint];
            target_s = (double) map_waypoints_s[target_waypoint];
          }
          double target_d = sqrtf(target_dx*target_dx + target_dy*target_dy);

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          // search around the road in the next 5 seconds to find if a lane change is necessary

          bool change_lane = false;
          // first check if there is a car in front in the same lane and whether a collision might occur
          for (int i = 0; i < sensor_fusion.size(); i++) {
            double obj_d = sensor_fusion[i][5];
            double obj_s = sensor_fusion[i][6];
            // check the car is in front and in the car's lane
            if (fabs(obj_d - car_d) < 1.2 && obj_s > car_s) {
              // check if maintaining current velocity will cause a collision
              double obj_yaw = sensor_fusion[i][4];
              double obj_v = sensor_fusion[i][3];
              double obj_v_s = sin(obj_yaw) * obj_v;
              double car_v_s = sin(car_yaw) * car_speed;
              if ((obj_v_s * t_end) + obj_s <= (car_v_s * t_end) + car_s) {
                // flag for lane change
                change_lane = true;
              }
              else {
                // set straight line waypoints
              }
            }
          }

          if (change_lane == true) {
            for (int i = 0; i < sensor_fusion.size(); i++) {
              double obj_d = (double) sensor_fusion[i][5];
              double obj_s = (double) sensor_fusion[i][6];

            }
          }

          // solve jerk minimizing trajectory for proposed goal and start states
          vector<double> jmt_s = JMT({car_s, car_v_s, 0}, {target_s, target_s_v, 1}, t_end);
          vector<double> jmt_d = JMT({car_d, car_v_d, t_end}, {target_d, target_d_v, 0.1}, t_end);

          // interpolate path points
          double n_s; double n_d; double n_t; vector<double> XY;
          for (int i = 0; i < path_vals; i++) {
            n_t = 0.02 * i;
            n_s = jmt_s[0] + jmt_s[1] * n_t + jmt_s[2] * (n_t*n_t) + jmt_s[3] * (n_t*n_t*n_t) +
                jmt_s[4] * (n_t*n_t*n_t*n_t) + jmt_s[5] * (n_t*n_t*n_t*n_t*n_t);
            n_d = jmt_d[0] + jmt_d[1] * n_t + jmt_d[2] * (n_t*n_t) + jmt_d[3] * (n_t*n_t*n_t) +
                  jmt_d[4] * (n_t*n_t*n_t*n_t) + jmt_d[5] * (n_t*n_t*n_t*n_t*n_t);
            XY = getXY(n_s, n_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            next_x_vals[i] = XY[0];
            next_y_vals[i] = XY[1];
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
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
















































































