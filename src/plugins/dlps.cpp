/* Copyright (c) 2017-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/host.h"
#include "simgrid/plugins/dlps.h"
#include "simgrid/plugins/dlps.hpp"
#include "simgrid/s4u/Link.hpp"
#include "src/surf/network_interface.hpp"
#include "src/surf/surf_interface.hpp"
#include "surf/surf.hpp"
#include "simgrid/sg_config.hpp"

#include <limits>

SIMGRID_REGISTER_PLUGIN(dlps, "Link DLPS.", &sg_dlps_plugin_init)

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(dlps, surf, "Logging specific to the SURF DLPS plugin");

namespace simgrid {
namespace plugin {

xbt::Extension<s4u::Link, DLPS> DLPS::EXTENSION_ID;

DLPS::DLPS(simgrid::s4u::Link* ptr) : link_(ptr), is_enabled_(false), dlps_mode_(simgrid::config::get_value<std::string>("network/dlps"))
{
  idle_threshold_laser_ = dlps_idle_threshold_laser;
  idle_threshold_tuning_ = dlps_idle_threshold_tuning;
  XBT_DEBUG("Instantiating a DLPS for link '%s'", link_->get_cname());
}

void DLPS::enable()
{
  xbt_assert(!is_enabled_, "Trying to enable load of link '%s' while it is already enabled, aborting.",
             link_->get_cname());
  XBT_DEBUG("Tracking load of link '%s'", link_->get_cname());

  is_enabled_ = true;
  // reset();
}

void DLPS::disable()
{
  xbt_assert(is_enabled_, "Trying to disable load of link '%s' while it is not enabled, aborting.", link_->get_cname());
  XBT_DEBUG("Untracking load of link '%s'", link_->get_cname());

  is_enabled_ = false;
}

void DLPS::reset()
{
  XBT_DEBUG("Resetting load of link '%s'", link_->get_cname());

  cumulated_bytes_      = 0.0;
  cumulated_energy_     = 0.0;
  min_bytes_per_second_ = std::numeric_limits<double>::max();
  max_bytes_per_second_ = std::numeric_limits<double>::lowest();
  XBT_DEBUG("min_bytes_per_second_ = %g", min_bytes_per_second_);
  XBT_DEBUG("max_bytes_per_second_ = %g", max_bytes_per_second_);
  last_reset_   = surf_get_clock();
  last_updated_ = last_reset_;
  idle_threshold_laser_ = 0.0;
  idle_threshold_tuning_ = 0.0;
  comm_trace.clear();
}

double DLPS::data_rate_to_power(double data_rate, bool laser_on, bool tuning_on)
{
  xbt_assert((data_rate > 0 && laser_on && tuning_on) || (data_rate == 0 && !(laser_on && !tuning_on)), "Something wrong!");
  double tuning_power = 0.96; // W
  double data_rate_dependent_power = data_rate / link_->get_bandwidth() * (0.63-0.03) + 0.03;

  return laser_on ? tuning_power + data_rate_dependent_power : (
         tuning_on ? tuning_power : 0.0);

}

void DLPS::update_on_comm_start(double actual_start_time)
{
  XBT_DEBUG("Updating load of link '%s' on communication start", link_->get_cname());
  xbt_assert(is_enabled_,
             "Trying to update load of link '%s' while it is NOT enabled, aborting."
             " Please enable your link with sg_dlps_enable before trying to access any of its load metrics.",
             link_->get_cname());

  std::string link_name = link_->get_cname();
  double current_instantaneous_bytes_per_second = link_->get_usage();
  double now                                    = surf_get_clock();

  s4u::Link::State last_state = link_->get_last_state();

  // Update minimum/maximum observed values if needed
  min_bytes_per_second_ = std::min(min_bytes_per_second_, current_instantaneous_bytes_per_second);
  max_bytes_per_second_ = std::max(max_bytes_per_second_, current_instantaneous_bytes_per_second);

  // Update cumulated load
  double duration_since_last_update = now - last_updated_;
  double bytes_since_last_update    = 0.;
  XBT_DEBUG("Cumulated %g bytes since last update (duration of %g seconds)", bytes_since_last_update,
            duration_since_last_update);
  cumulated_bytes_ += bytes_since_last_update;

  // Update cumulated energy
  double current_instantaneous_power = 0.0;
  double energy_since_last_update    = 0.0;

  if (last_updated_ < 0) { // if it is the first communication, "on-off" and "full" will be OFF before this
    xbt_assert(current_instantaneous_bytes_per_second == 0, "Link usage should be 0 if this is the first communication.");
    if (dlps_mode_ == "none") {
      // energy
      current_instantaneous_power = data_rate_to_power(link_->get_bandwidth() * sg_bandwidth_factor);
      energy_since_last_update    = now * current_instantaneous_power;
      // trace
//      comm_trace.push_back(std::make_tuple(0.0, "ON", 0.0, 0.0, 0.0));
//      comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, current_instantaneous_power, cumulated_energy_ + energy_since_last_update));
    } else if (dlps_mode_ == "laser") {
      // energy
      current_instantaneous_power = data_rate_to_power(0.0, false, true);
      energy_since_last_update    = now * current_instantaneous_power;
      // trace
//      comm_trace.push_back(std::make_tuple(0.0, "STANDBY", 0.0, 0.0, 0.0));
//      comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, current_instantaneous_power, cumulated_energy_ + energy_since_last_update));
    } else {
      // no energy consumption before the first communication
      // trace
//      comm_trace.push_back(std::make_tuple(0.0, "OFF", 0.0, 0.0, 0.0));
//      comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, 0.0, 0.0));
    }
  } else { // last_updated_ >= 0
    if (link_->get_num_active_actions_at(now) == 1 && now > last_updated_) { // The first of several consecutive start, compute the energy after last end
      xbt_assert(current_instantaneous_bytes_per_second == 0, "Link usage should be 0 if this is the first of several consecutive communication.");
      if (dlps_mode_ == "full") {
	// ready energy
        double ready_time = std::min(duration_since_last_update, idle_threshold_laser_);
        double ready_power = data_rate_to_power(0.0, true, true);
        energy_since_last_update += ready_time * ready_power;
        // trace
	if (duration_since_last_update <= idle_threshold_laser_) { // idle not long enough to go into standby
//	  comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, ready_power, cumulated_energy_ + energy_since_last_update));
	} else { // idle long enough to go into standby
//	  comm_trace.push_back(std::make_tuple(last_updated_ + idle_threshold_laser_, "STANDBY", 0.0, ready_power, cumulated_energy_ + energy_since_last_update));
	  // standby energy
          double standby_time = std::min(duration_since_last_update - idle_threshold_laser_, idle_threshold_tuning_ - idle_threshold_laser_);
	  double standby_power = data_rate_to_power(0.0, false, true);
          energy_since_last_update += standby_time * standby_power;
	  // trace
	  if (duration_since_last_update <= idle_threshold_tuning_) { // idle not long enough to go into off
//	    comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, standby_power, cumulated_energy_ + energy_since_last_update));
	  } else { // idle long enough to go into off
//	    comm_trace.push_back(std::make_tuple(last_updated_ + idle_threshold_tuning_, "OFF", 0.0, standby_power, cumulated_energy_ + energy_since_last_update));
//	    comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, 0.0, cumulated_energy_ + energy_since_last_update));
	  }
	}
      } else if (dlps_mode_ == "laser") {
	// ready energy
        double ready_time = std::min(duration_since_last_update, idle_threshold_laser_);
	double ready_power = data_rate_to_power(0.0, true, true);
	energy_since_last_update += ready_time * ready_power;
	// trace
	if (duration_since_last_update <= idle_threshold_laser_) { // idle not long enough to go into standby
//	  comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, ready_power, cumulated_energy_ + energy_since_last_update));
	} else { // idle long enough to go into standby
//	  comm_trace.push_back(std::make_tuple(last_updated_ + idle_threshold_laser_, "STANDBY", 0.0, ready_power, cumulated_energy_ + energy_since_last_update));
	  // standby energy
          double standby_time = std::max(duration_since_last_update - idle_threshold_laser_, 0.0);
	  double standby_power = data_rate_to_power(0.0, false, true);
	  energy_since_last_update += standby_time * standby_power;
	  // trace
//	  comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, 0.0, cumulated_energy_ + energy_since_last_update));
	}
      } else if (dlps_mode_ == "on-off") {
        // no energy after last end
//	comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, 0.0, cumulated_energy_));
      } else {
        current_instantaneous_power = data_rate_to_power(link_->get_bandwidth() * sg_bandwidth_factor);
        energy_since_last_update    = duration_since_last_update * current_instantaneous_power;
//	comm_trace.push_back(std::make_tuple(now, "STARTING", 0.0, current_instantaneous_power, cumulated_energy_ + energy_since_last_update));
      }
    } else { // Not the first of several consecutive start
      current_instantaneous_power = data_rate_to_power(
        (dlps_mode_ == "full" || dlps_mode_ == "laser" || dlps_mode_ == "on-off") ?
        current_instantaneous_bytes_per_second : link_->get_bandwidth() * sg_bandwidth_factor);
      energy_since_last_update    = duration_since_last_update * current_instantaneous_power;
//      comm_trace.push_back(std::make_tuple(now, "STARTING", current_instantaneous_bytes_per_second, current_instantaneous_power, cumulated_energy_ + energy_since_last_update));
    }
  }

  comm_trace.push_back(std::make_tuple(link_->get_next_wake(), link_->get_num_active_actions_at(now), current_instantaneous_bytes_per_second, current_instantaneous_power));

  XBT_DEBUG("Cumulated %g J since last update (duration of %g seconds)", energy_since_last_update,
            duration_since_last_update);
  cumulated_energy_ += energy_since_last_update;

  // Update idle thresholds based on values in the circular buffer
  /* if (dlps_mode_ == "full" || dlps_mode_ == "laser") {
    
    idle_threshold_laser_ = link_->interval_recorder[0] > dlps_idle_threshold_laser && link_->interval_recorder[1] > dlps_idle_threshold_laser ? 0.0 : (
                             link_->interval_recorder[0] > idle_threshold_laser_ && link_->interval_recorder[1] > idle_threshold_laser_ ? link_->interval_recorder[1] : (
                             link_->interval_recorder[0] < idle_threshold_laser_ && link_->interval_recorder[1] < idle_threshold_laser_ ? link_->interval_recorder[1] : idle_threshold_laser_));

    if (dlps_mode_ == "full") {
      idle_threshold_tuning_ = link_->interval_recorder[0] > dlps_idle_threshold_tuning && link_->interval_recorder[1] > dlps_idle_threshold_tuning ? idle_threshold_laser_ : (
                               link_->interval_recorder[0] > idle_threshold_tuning_ && link_->interval_recorder[1] > idle_threshold_tuning_ ? link_->interval_recorder[1] : (
                               link_->interval_recorder[0] < idle_threshold_tuning_ && link_->interval_recorder[1] < idle_threshold_tuning_ ? std::max(link_->interval_recorder[1], idle_threshold_laser_) : idle_threshold_tuning_));

      if (idle_threshold_tuning_ < idle_threshold_laser_)
	idle_threshold_tuning_ = idle_threshold_laser_;
    }
  } */

  last_updated_ = now;

  xbt_assert(bytes_since_last_update >= 0, "DLPS plugin inconsistency: negative amount of bytes is accumulated.");
}

void DLPS::update_on_comm_end(double actual_start_time, double actual_transfer_end_time, double size)
{
  XBT_DEBUG("Updating link '%s' on communication end", link_->get_cname());
  xbt_assert(is_enabled_,
             "Trying to update link '%s' while DLPS is NOT enabled, aborting."
             " Please enable your link with sg_dlps_enable before trying to access any of its metrics.",
             link_->get_cname());

  std::string link_name = link_->get_cname();
  double current_instantaneous_bytes_per_second = link_->get_usage();
  double now                                    = actual_transfer_end_time;

  s4u::Link::State last_state = link_->get_last_state();

  switch (last_state) {
    case s4u::Link::State::OFF:
      xbt_assert(dlps_mode_ == "full" || dlps_mode_ == "on-off", "Link state is OFF which cannot happen for current DLPS mode.");
      break;
    case s4u::Link::State::STANDBY:
      xbt_assert(dlps_mode_ == "full" || dlps_mode_ == "laser", "Link state is STANDBY which cannot happen for current DLPS mode.");
      break;
    case s4u::Link::State::READY:
      xbt_assert(dlps_mode_ == "full" || dlps_mode_ == "laser", "Link state is READY which cannot happen for current DLPS mode.");
      break;
  }

  // Update minimum/maximum observed values if needed
  min_bytes_per_second_ = std::min(min_bytes_per_second_, current_instantaneous_bytes_per_second);
  max_bytes_per_second_ = std::max(max_bytes_per_second_, current_instantaneous_bytes_per_second);

  // Update cumulated load
  double duration_since_last_update = now - std::max(last_updated_, actual_start_time);
  double bytes_since_last_update    = size;
  XBT_DEBUG("Cumulated %g bytes since last update (duration of %g seconds)", bytes_since_last_update,
            duration_since_last_update);
  cumulated_bytes_ += bytes_since_last_update;

  // Update cumulated energy
  double current_instantaneous_power = 0.0;
  double energy_since_last_update    = 0.0;

  current_instantaneous_power = data_rate_to_power(
      (dlps_mode_ == "full" || dlps_mode_ == "laser" || dlps_mode_ == "on-off") ?
      current_instantaneous_bytes_per_second : link_->get_bandwidth() * sg_bandwidth_factor);

  // delay energy
  double delay_time = std::max(actual_start_time - last_updated_, 0.0);
  energy_since_last_update += current_instantaneous_power * delay_time;
  // trace
//  comm_trace.push_back(std::make_tuple(last_updated_ + delay_time, "STARTED", current_instantaneous_bytes_per_second, current_instantaneous_power, cumulated_energy_ + energy_since_last_update));
  // actual transfer energy
  double transfer_time = now - actual_start_time;
  energy_since_last_update += current_instantaneous_power * transfer_time;
  // trace
//  std::string to_state = (dlps_mode_ == "full" || dlps_mode_ == "laser") ? "READY" : (
//			  dlps_mode_ == "on-off" ? "OFF" : "ON");
//  comm_trace.push_back(std::make_tuple(now, to_state, current_instantaneous_bytes_per_second, current_instantaneous_power, cumulated_energy_ + energy_since_last_update));

  comm_trace.push_back(std::make_tuple(actual_transfer_end_time, link_->get_num_active_actions_at(actual_transfer_end_time), current_instantaneous_bytes_per_second, current_instantaneous_power));


  XBT_DEBUG("Cumulated %g J since last update (duration of %g seconds)", energy_since_last_update,
            duration_since_last_update);
  cumulated_energy_ += energy_since_last_update;

  last_updated_ = now;

  xbt_assert(bytes_since_last_update >= 0, "DLPS plugin inconsistency: negative amount of bytes is accumulated.");
}

s4u::Link* DLPS::get_s4u_link() {
  return link_;
}

bool DLPS::is_enabled() const
{
  return is_enabled_;
}

std::string DLPS::get_dlps_mode() const
{
  return dlps_mode_;
}

double DLPS::get_last_updated()
{
  return last_updated_;
}

double DLPS::get_idle_threshold_laser()
{
  return idle_threshold_laser_;
}

double DLPS::get_idle_threshold_tuning()
{
  return idle_threshold_tuning_;
}

double DLPS::get_cumulated_bytes()
{
  return cumulated_bytes_;
}

double DLPS::get_cumulated_energy()
{
  update_on_comm_start(surf_get_clock());
  return cumulated_energy_;
}

double DLPS::get_min_bytes_per_second()
{
  return min_bytes_per_second_;
}

double DLPS::get_max_bytes_per_second()
{
  return max_bytes_per_second_;
}

double DLPS::get_average_bytes()
{
  double now = surf_get_clock();
  if (now > last_reset_)
    return cumulated_bytes_ / (now - last_reset_);
  else
    return 0;
}

} // namespace plugin
} // namespace simgrid

using simgrid::plugin::DLPS;

/* **************************** events  callback *************************** */
static void on_communicate(simgrid::kernel::resource::NetworkAction& action)
{
  // XBT_DEBUG("on_communicate is called");
  
  double now = surf_get_clock();
  for (auto* link : action.get_route()) {
    if (link != nullptr && link->get_sharing_policy() != simgrid::s4u::Link::SharingPolicy::WIFI) {
      auto dlps = link->get_iface()->extension<DLPS>();
      if (dlps->is_enabled()) {
        XBT_INFO("%.17f,%ld\n", now, link->get_iface()->get_num_active_actions_at(now));
        dlps->update_on_comm_start(action.get_actual_start_time());
      }
    }
  }
}

static void on_communication_state_change(const simgrid::kernel::resource::NetworkAction& action,
         simgrid::kernel::resource::Action::State previous_state) {

  double now = surf_get_clock();
  double transfer_time_per_link = (now - action.get_actual_start_time()) / action.get_route().size();
  int link_idx = 0;

  for (auto* link : action.get_route()) {
    link_idx++;
    if (link != nullptr && link->get_sharing_policy() != simgrid::s4u::Link::SharingPolicy::WIFI) {
      auto dlps = link->get_iface()->extension<DLPS>();
      if (dlps->is_enabled()) {
        if (action.get_state() == simgrid::kernel::resource::Action::State::FINISHED) {

	  double actual_transfer_end_time = action.get_actual_start_time() + link_idx * transfer_time_per_link;
          link->get_iface()->remove_active_action_at(action.get_start_time());
          link->get_iface()->set_last_busy(actual_transfer_end_time);

          if (link->get_iface()->get_num_active_actions_at(actual_transfer_end_time) == 0) {
            if (dlps->get_dlps_mode() == "full") {
              link->get_iface()->set_next_ready(actual_transfer_end_time);
              link->get_iface()->set_next_standby(actual_transfer_end_time + dlps->get_idle_threshold_laser());
              link->get_iface()->set_next_off(actual_transfer_end_time + dlps->get_idle_threshold_tuning());
            }
            else if (dlps->get_dlps_mode() == "laser") {
              link->get_iface()->set_next_ready(actual_transfer_end_time);
              link->get_iface()->set_next_standby(actual_transfer_end_time + dlps->get_idle_threshold_tuning());
            } else if (dlps->get_dlps_mode() == "on-off") {
              link->get_iface()->set_next_off(actual_transfer_end_time);
            }
          }
          dlps->update_on_comm_end(action.get_actual_start_time(), actual_transfer_end_time, action.get_size());
	  XBT_INFO("%.17f,%ld\n", now, link->get_iface()->get_num_active_actions_at(now));
	}
      }
    }
  }
}

/* **************************** Public interface *************************** */

/**
 * @ingroup plugin_dlps
 * @brief Initialize the link cumulated load plugin.
 * @pre The energy plugin should NOT be initialized.
 */
void sg_dlps_plugin_init()
{
  if (DLPS::EXTENSION_ID.valid()) {
    return;
  }

  DLPS::EXTENSION_ID = simgrid::s4u::Link::extension_create<DLPS>();

  if (simgrid::s4u::Engine::is_initialized()) {
    const simgrid::s4u::Engine* e = simgrid::s4u::Engine::get_instance();
    for (auto& link : e->get_all_links()) {
      if (link->get_sharing_policy() != simgrid::s4u::Link::SharingPolicy::WIFI) {
        link->extension_set(new DLPS(link));
      }
    }
  }

  // Attach new DLPS to links created in the future.
  simgrid::s4u::Link::on_creation.connect([](simgrid::s4u::Link& link) {
    if (link.get_sharing_policy() != simgrid::s4u::Link::SharingPolicy::WIFI) {
      XBT_DEBUG("Wired link '%s' created. Attaching DLPS to it.", link.get_cname());
      link.extension_set(new DLPS(&link));
    } else {
      XBT_DEBUG("Wireless link '%s' created. NOT attaching any DLPS to it.", link.get_cname());
    }
  });

  // Call this plugin on some of the links' events.
  simgrid::s4u::Link::on_communication_state_change.connect(&on_communication_state_change);
  simgrid::s4u::Link::on_communicate.connect(&on_communicate);
}

/**
 * @ingroup plugin_dlps
 * @brief Start the attaching of a link.
 * @details This is required so the link cumulated load can be obtained later on.
 * @pre The energy plugin should be initialized (cf. sg_dlps_plugin_init()).
 * @pre The link should be in "untracked" state. In other words, do not call this function twice on the same link
 * without a sg_dlps_untrack() call between them.
 *
 * @param link The link to enable.
 */
void sg_dlps_enable(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_enable. Aborting.");
  link->extension<DLPS>()->enable();
}

/**
 * @ingroup plugin_dlps
 * @brief Stop the tracking of a link.
 * @details Once the tracking is stopped, the cumulated load of the link can no longer be obtained until
 * sg_dlps_track() is called again on this link.
 * @pre The energy plugin should be initialized (cf. sg_dlps_plugin_init()).
 * @pre The link should be in "tracked" state. In other words, do not call this function twice on the same link without
 * a sg_dlps_track() call between them.
 *
 * @param link The link to untrack.
 */
void sg_dlps_disable(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_disable. Aborting.");
  link->extension<DLPS>()->disable();
}

/**
 * @ingroup plugin_dlps
 * @brief Resets the cumulated load counters of a link.
 * @pre The energy plugin should be initialized (cf. sg_dlps_plugin_init()).
 * @pre The link should be in "tracked" state (cf. sg_dlps_track()).
 *
 * @param link The link whose counters should be reset.
 */
void sg_dlps_reset(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_reset. Aborting.");
  link->extension<DLPS>()->reset();
}

/**
 * @ingroup plugin_dlps
 * @brief Get the cumulated load of a link (since the last call to sg_dlps_reset()).
 * @pre The energy plugin should be initialized (cf. sg_dlps_plugin_init()).
 * @pre The link should be in "tracked" state (cf. sg_dlps_track()).

 * @param link The link whose cumulated load is requested.
 * @return The load (in bytes) that passed through the given link since the last call to sg_dlps_reset.
 */
double sg_dlps_get_cum_load(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_get_cum_load. Aborting.");
  return link->extension<DLPS>()->get_cumulated_bytes();
}

double sg_dlps_get_cum_energy(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_get_cum_load. Aborting.");
  return link->extension<DLPS>()->get_cumulated_energy();
}

/**
 * @ingroup plugin_dlps
 * @brief Get the average load of a link (since the last call to sg_dlps_reset()).
 * @pre The energy plugin should be initialized (cf. sg_dlps_plugin_init()).
 * @pre The link should be in "tracked" state (cf. sg_dlps_track()).

 * @param link The link whose average load is requested.
 * @return The average load (in bytes) that passed through the given link since the last call to sg_dlps_reset.
 */
double sg_dlps_get_avg_load(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_get_avg_load. Aborting.");
  return link->extension<DLPS>()->get_average_bytes();
}

/**
 * @ingroup plugin_dlps
 * @brief Get the minimum instantaneous load of a link (since the last call to sg_dlps_reset()).
 * @pre The energy plugin should be initialized (cf. sg_dlps_plugin_init()).
 * @pre The link should be in "tracked" state (cf. sg_dlps_track()).

 * @param link The link whose average load is requested.
 * @return The minimum load instantaneous load (in bytes per second) that passed through the given link since the last
 call to sg_dlps_reset.
 */
double sg_dlps_get_min_instantaneous_load(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_get_min_load. Aborting.");
  return link->extension<DLPS>()->get_min_bytes_per_second();
}

/**
 * @ingroup plugin_dlps
 * @brief Get the maximum instantaneous load of a link (since the last call to sg_dlps_reset()).
 * @pre The energy plugin should be initialized (cf. sg_dlps_plugin_init()).
 * @pre The link should be in "tracked" state (cf. sg_dlps_track()).

 * @param link The link whose average load is requested.
 * @return The maximum load instantaneous load (in bytes per second) that passed through the given link since the last
 call to sg_dlps_reset.
 */
double sg_dlps_get_max_instantaneous_load(const_sg_link_t link)
{
  xbt_assert(DLPS::EXTENSION_ID.valid(),
             "Please call sg_dlps_plugin_init before sg_dlps_get_max_load. Aborting.");
  return link->extension<DLPS>()->get_max_bytes_per_second();
}
