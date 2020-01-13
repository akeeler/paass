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

  std::vector<std::pair<double, double>> gammaEvents = {};
  std::vector<std::pair<double, double>> gamma_blue = {};
  std::vector<std::pair<double, double>> gamma_black = {};
  std::vector<std::pair<double, double>> gamma_green = {};
  std::vector<std::pair<double, double>> gamma_red = {};

};
#endif
