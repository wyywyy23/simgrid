/* Copyright (c) 2004-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef S4U_LINK_HPP
#define S4U_LINK_HPP

#include <simgrid/forward.h>
#include <simgrid/kernel/resource/Action.hpp>
#include <simgrid/link.h>
#include <string>
#include <unordered_map>
#include <xbt/Extendable.hpp>
#include <xbt/base.h>
#include <xbt/signal.hpp>
#include <boost/circular_buffer.hpp>

/***********
 * Classes *
 ***********/

namespace simgrid {

extern template class XBT_PUBLIC xbt::Extendable<s4u::Link>;

namespace s4u {
/**
 * @beginrst
 * A Link represents the network facilities between :cpp:class:`hosts <simgrid::s4u::Host>`.
 * @endrst
 */
class XBT_PUBLIC Link : public xbt::Extendable<Link> {
  friend kernel::resource::LinkImpl;

  // Links are created from the NetZone, and destroyed by their private implementation when the simulation ends
  explicit Link(kernel::resource::LinkImpl* pimpl) : pimpl_(pimpl), interval_recorder(2) {}
  virtual ~Link() = default;
  // The private implementation, that never changes
  kernel::resource::LinkImpl* const pimpl_;

public:
  enum class SharingPolicy { WIFI = 3, SPLITDUPLEX = 2, SHARED = 1, FATPIPE = 0 };

  kernel::resource::LinkImpl* get_impl() const { return pimpl_; }

  /** @brief Retrieve a link from its name */
  static Link* by_name(const std::string& name);
  static Link* by_name_or_null(const std::string& name);

  /** @brief Retrieves the name of that link as a C++ string */
  const std::string& get_name() const;
  /** @brief Retrieves the name of that link as a C string */
  const char* get_cname() const;

  /** Get/Set the bandwidth of the current Link (in bytes per second) */
  double get_bandwidth() const;
  Link* set_bandwidth(double value);

  /** Get/Set the latency of the current Link (in seconds) */
  double get_latency() const;
  /**
   * @brief Set link's latency
   *
   * @param value New latency value (in s)
   */
  Link* set_latency(double value);
  /**
   * @brief Set latency (string version)
   *
   * @throw std::invalid_argument if latency format is incorrect.
   */
  Link* set_latency(const std::string& value);

  /** @brief Describes how the link is shared between flows */
  Link* set_sharing_policy(SharingPolicy policy);
  SharingPolicy get_sharing_policy() const;

  /** Setup the profile with states events (ON or OFF). The profile must contain boolean values. */
  Link* set_state_profile(kernel::profile::Profile* profile);
  /** Setup the profile with bandwidth events (peak speed changes due to external load).
   * The profile must contain percentages (value between 0 and 1). */
  Link* set_bandwidth_profile(kernel::profile::Profile* profile);
  /** Setup the profile file with latency events (peak latency changes due to external load).
   * The profile must contain absolute values */
  Link* set_latency_profile(kernel::profile::Profile* profile);

  const std::unordered_map<std::string, std::string>* get_properties() const;
  const char* get_property(const std::string& key) const;
  Link* set_properties(const std::unordered_map<std::string, std::string>& properties);
  Link* set_property(const std::string& key, const std::string& value);

  /** @brief Set the level of communication speed of the given host on this wifi link.
   *
   * The bandwidth of a wifi link for a given host depends on its SNR (signal to noise ratio),
   * which ultimately depends on the distance between the host and the station and the material between them.
   *
   * This is modeled in SimGrid by providing several bandwidths to wifi links, one per SNR level (just provide
   * comma-separated values in the XML file). By default, the first level in the list is used, but you can use the
   * current function to specify that a given host uses another level of bandwidth. This can be used to take the
   * location of hosts into account, or even to model mobility in your SimGrid simulation.
   *
   * Note that this function asserts that the link is actually a wifi link */
  void set_host_wifi_rate(const s4u::Host* host, int level) const;

  /** @brief Returns the current load (in bytes per second) */
  double get_usage() const;

  /** @brief Check if the Link is used (at least one flow uses the link) */
  bool is_used() const;

  /** @brief Check if the Link is shared (not a FATPIPE) */
  bool is_shared() const;

  void turn_on();
  void turn_off();

  /** wyy: more states for links */
  XBT_DECLARE_ENUM_CLASS(State, OFF, STANDBY, READY, ON);
  void set_last_state(Link::State state) { last_state_ = state; }
  Link::State get_last_state() const { return last_state_; }
  const char* get_last_state_str() const;
  Link::State last_state_ = Link::State::ON;

  /** wyy: timestamps related to the states */
  double next_wake_ = 0.0;
  void set_next_wake(double time) { next_wake_ = time; }
  double get_next_wake() const { return next_wake_; }

  double next_on_ = 0.0;
  void set_next_on(double time) { next_on_ = time; }
  double get_next_on() const { return next_on_; }

  double next_ready_ = std::numeric_limits<double>::max();
  void set_next_ready(double time) { next_ready_ = time; }
  double get_next_ready() const { return next_ready_; }

  double next_standby_ = std::numeric_limits<double>::max();
  void set_next_standby(double time) { next_standby_ = time; }
  double get_next_standby() const { return next_standby_; }

  double next_off_ = std::numeric_limits<double>::max();
  void set_next_off(double time) { next_off_ = time; }
  double get_next_off() const { return next_off_; }

  double last_busy_ = -1.0;
  void set_last_busy(double time) { last_busy_ = time; }
  double get_last_busy() const { return last_busy_; }

  boost::circular_buffer<double> interval_recorder;

  /** Keep track of active actions using this link */
  std::map<double, unsigned long> active_actions;
  unsigned long get_num_active_actions_before(double time);
  unsigned long get_num_active_actions_at(double time);
  void add_active_action_at(double time);
  void remove_active_action_at(double time);

  bool is_on() const;

  void seal();

  /* The signals */
  /** @brief Callback signal fired when a new Link is created */
  static xbt::signal<void(Link&)> on_creation;

  /** @brief Callback signal fired when a Link is destroyed */
  static xbt::signal<void(Link const&)> on_destruction;

  /** @brief Callback signal fired when the state of a Link changes (when it is turned on or off) */
  static xbt::signal<void(Link const&)> on_state_change;

  /** @brief Callback signal fired when the bandwidth of a Link changes */
  static xbt::signal<void(Link const&)> on_bandwidth_change;

  /** @brief Callback signal fired when a communication starts */
  static xbt::signal<void(kernel::resource::NetworkAction&)> on_communicate;

  /** @brief Callback signal fired when a communication changes it state (ready/done/cancel) */
  static xbt::signal<void(kernel::resource::NetworkAction&, kernel::resource::Action::State)>
      on_communication_state_change;

  /** Unified signal for DLPS <link, last_state, total_rate>*/
  static xbt::signal<void(Link&, Link::State, double)> on_dlps_state_change;
};
} // namespace s4u
} // namespace simgrid

#endif /* S4U_LINK_HPP */
