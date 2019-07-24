#ifndef __PIDSTRUCT_HPP__
#define __PIDSTRUCT_HPP__

struct PidEvent {

  bool pid = false;
  double pin = -1;
  double beam_tof = -1;
  double i2pos = -10000;
  double corrected_beam_tof = -1;
};
#endif
