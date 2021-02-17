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
  void update_on_comm_start();
  void update_on_comm_end();
  void update_load();

  /// Getter methods.
  bool is_enabled() const;
  s4u::Link* get_s4u_link();
  double get_last_updated();
  double get_average_bytes();
  double get_cumulated_bytes();
  double get_cumulated_energy();
  double get_min_bytes_per_second();
  double get_max_bytes_per_second();

private:
  s4u::Link* link_{};      /*< The link onto which this data is enabled*/
  bool is_enabled_{false}; /*<Whether the link is enabled or not*/

  double original_latency_{};  // Back up of original latency;

  double cumulated_bytes_{0.0};      /*< Cumulated load since last reset*/
  double min_bytes_per_second_{0.0}; /*< Minimum instantaneous load observed since last reset*/
  double max_bytes_per_second_{0.0}; /*< Maximum instantaneous load observed since last reset*/
  double last_reset_{-1.0};          /*< Timestamp of the last reset (init timestamp by default)*/
  double last_updated_{-1.0};        /*< Timestamp of the last update event*/
};

}
}

#endif
