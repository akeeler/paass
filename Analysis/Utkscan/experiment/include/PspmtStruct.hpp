#ifndef __PSPMTSTRUCT_HPP__
#define __PSPMTSTRUCT_HPP__

#include "PidStruct.hpp"

struct PspmtEvent {

  PidEvent pidEvent = {};

  double x_position = -2;
  double y_position = -2;
  int x_pixel = -2;
  int y_pixel = -2;
  int pixel_num = -2;
  double event_time = -2;
  bool implant = false;
  bool decay = false;
  std::pair<double, double> low_dynode = {};
  std::pair<double, double> hi_dynode = {};
  double low_dynode_mult = 0;
  double hi_dynode_mult = 0;
  std::vector<std::pair<double, double>> gammaEvents = {};

};
#endif
