/* Copyright (c) 2017-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SIMGRID_PLUGINS_DLPS_HPP_
#define SIMGRID_PLUGINS_DLPS_HPP_

#include "simgrid/s4u.hpp"

namespace simgrid {
namespace plugin {

class DLPS {
public:
  static xbt::Extension<s4u::Link, DLPS> EXTENSION_ID;

  explicit DLPS(simgrid::s4u::Link* ptr);
  ~DLPS() = default;

  void enable();
  void disable();
  void reset();
  void update_on_comm_start(unsigned long, double);
  void update_on_comm_end(unsigned long, double, double);
  void update_load();

  /// Getter methods.
  bool is_enabled() const;
  std::string get_dlps_mode() const;
  s4u::Link* get_s4u_link();
  double get_last_updated();
  double get_idle_threshold_laser();
  double get_idle_threshold_tuning();
  double get_average_bytes();
  double get_cumulated_bytes();
  double get_cumulated_energy();
  double get_min_bytes_per_second();
  double get_max_bytes_per_second();

  std::vector<std::tuple<double, unsigned long, std::string, double, double, double>> get_comm_trace() const { return comm_trace; }

private:
  s4u::Link* link_{};      /*< The link onto which this data is enabled*/
  bool is_enabled_{false}; /*< Whether the link is enabled or not*/
  std::string dlps_mode_{"none"}; /*< DLPS mode */

  double cumulated_bytes_{0.0};      /*< Cumulated load since last reset*/
  double cumulated_energy_{0.0};     /*< Cumulated energy since last reset*/
  double min_bytes_per_second_{std::numeric_limits<double>::max()}; /*< Minimum instantaneous load observed since last reset*/
  double max_bytes_per_second_{std::numeric_limits<double>::lowest()}; /*< Maximum instantaneous load observed since last reset*/
  double last_reset_{0.0};          /*< Timestamp of the last reset (init timestamp by default)*/
  double last_updated_{-1.0};        /*< Timestamp of the last update event*/

  double idle_threshold_laser_{0.0};  /*< Idle threshold laser for the managed link*/
  double idle_threshold_tuning_{0.0}; /*< Idle threshold tuning for the managed link*/

  double data_rate_to_power(double rate, bool laser_on = true, bool tuning_on = true); /*< Compute power from data rate*/

  std::vector<std::tuple<double, unsigned long, std::string, double, double, double>> comm_trace; /*< now, action_id, state, action_actual_start/finish, link_catering_start/finish, rate */

};

}
}

#endif
