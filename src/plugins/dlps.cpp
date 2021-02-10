/* Copyright (c) 2017-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/host.h"
#include "simgrid/plugins/dlps.h"
#include "simgrid/s4u/Link.hpp"
#include "src/surf/network_interface.hpp"
#include "src/surf/surf_interface.hpp"
#include "surf/surf.hpp"
#include "simgrid/s4u.hpp"

#include <limits>

SIMGRID_REGISTER_PLUGIN(dlps, "Link DLPS.", &sg_dlps_plugin_init)

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(dlps, surf, "Logging specific to the SURF DLPS plugin");

const double delay_tuning = 1.0e-3;
const double delay_laser_stabilizing = 10.0e-9;
const double delay_laser_waking = 1.0e-9;

const double idle_threshold_tuning = 1.0e-3;
const double idle_threshold_laser = 300.0e-9;

namespace simgrid {
namespace plugin {

class DLPS {
public:
  static simgrid::xbt::Extension<simgrid::s4u::Link, DLPS> EXTENSION_ID;

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
  double get_last_busy();
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
  double last_busy_{-1.0};            // Timestamp when the last communication ended
};

xbt::Extension<s4u::Link, DLPS> DLPS::EXTENSION_ID;

DLPS::DLPS(simgrid::s4u::Link* ptr) : link_(ptr), is_enabled_(false), original_latency_(ptr->get_latency())
{
  XBT_DEBUG("Instantiating a DLPS for link '%s'", link_->get_cname());
}

void DLPS::enable()
{
  xbt_assert(!is_enabled_, "Trying to enable load of link '%s' while it is already enabled, aborting.",
             link_->get_cname());
  XBT_DEBUG("Tracking load of link '%s'", link_->get_cname());

  is_enabled_ = true;
  reset();
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
  min_bytes_per_second_ = std::numeric_limits<double>::max();
  max_bytes_per_second_ = std::numeric_limits<double>::lowest();
  XBT_DEBUG("min_bytes_per_second_ = %g", min_bytes_per_second_);
  XBT_DEBUG("max_bytes_per_second_ = %g", max_bytes_per_second_);
  last_reset_   = surf_get_clock();
  last_updated_ = last_reset_;
}

void DLPS::update_load()
{
  XBT_DEBUG("Updating load of link '%s'", link_->get_cname());
  xbt_assert(is_enabled_,
             "Trying to update load of link '%s' while it is NOT enabled, aborting."
             " Please enable your link with sg_dlps_enable before trying to access any of its load metrics.",
             link_->get_cname());

  std::string link_name = link_->get_cname();
  double current_instantaneous_bytes_per_second = link_->get_usage();
  double now                                    = surf_get_clock();

  // Update minimum/maximum observed values if needed
  min_bytes_per_second_ = std::min(min_bytes_per_second_, current_instantaneous_bytes_per_second);
  max_bytes_per_second_ = std::max(max_bytes_per_second_, current_instantaneous_bytes_per_second);

  // Update cumulated load
  double duration_since_last_update = now - last_updated_;
  double bytes_since_last_update    = duration_since_last_update * current_instantaneous_bytes_per_second;
  XBT_DEBUG("Cumulated %g bytes since last update (duration of %g seconds)", bytes_since_last_update,
            duration_since_last_update);
  xbt_assert(bytes_since_last_update >= 0, "DLPS plugin inconsistency: negative amount of bytes is accumulated.");

  cumulated_bytes_ += bytes_since_last_update;
  last_updated_ = now;
}

void DLPS::update_on_comm_start()
{
  XBT_DEBUG("Updating load of link '%s' on communication start", link_->get_cname());
  xbt_assert(is_enabled_,
             "Trying to update load of link '%s' while it is NOT enabled, aborting."
             " Please enable your link with sg_dlps_enable before trying to access any of its load metrics.",
             link_->get_cname());

  std::string link_name = link_->get_cname();
  double current_instantaneous_bytes_per_second = link_->get_usage();
  double now                                    = surf_get_clock();
  double this_latency = original_latency_;

  // Check the interval from last busy
  if (last_busy_ < 0) { // Never communicated
    this_latency += delay_tuning + delay_laser_stabilizing;
    link_->set_latency(this_latency);
    XBT_INFO("%s,tuning,%.17f,%f\n", link_name.c_str(), now, current_instantaneous_bytes_per_second);
    XBT_INFO("%s,lasering,%.17f,%f\n", link_name.c_str(), now + delay_tuning, current_instantaneous_bytes_per_second);
    XBT_INFO("%s,communicate,%.17f,%f\n", link_name.c_str(), now + delay_tuning + delay_laser_stabilizing, current_instantaneous_bytes_per_second);
  } else if (now - last_busy_ > idle_threshold_tuning) { // Idled till tuning is off
    this_latency += delay_tuning + delay_laser_stabilizing;
    link_->set_latency(this_latency);
    XBT_INFO("%s,tuning,%.17f,%f\n", link_name.c_str(), now, current_instantaneous_bytes_per_second);
    XBT_INFO("%s,lasering,%.17f,%f\n", link_name.c_str(), now + delay_tuning, current_instantaneous_bytes_per_second);
    XBT_INFO("%s,communicate,%.17f,%f\n", link_name.c_str(), now + delay_tuning + delay_laser_stabilizing, current_instantaneous_bytes_per_second);
  } else if (now - last_busy_ > idle_threshold_laser) { // Idled till laser is off, but tuning is on
    this_latency += delay_laser_stabilizing;
    link_->set_latency(this_latency);
    XBT_INFO("%s,lasering,%.17f,%f\n", link_name.c_str(), now, current_instantaneous_bytes_per_second);
    XBT_INFO("%s,communicate,%.17f,%f\n", link_name.c_str(), now + delay_laser_stabilizing, current_instantaneous_bytes_per_second);
  } else if (now - last_busy_ > 0) { // Idled till laser is sleeping
    this_latency += delay_laser_waking;
    link_->set_latency(this_latency);
    XBT_INFO("%s,waking,%.17f,%f\n", link_name.c_str(), now, current_instantaneous_bytes_per_second);
    XBT_INFO("%s,communicate,%.17f,%f\n", link_name.c_str(), now + delay_laser_waking, current_instantaneous_bytes_per_second);
  } else {  // Did not sleep at all
    link_->set_latency(this_latency);
    XBT_INFO("%s,communicate,%.17f,%f\n", link_name.c_str(), now, current_instantaneous_bytes_per_second);
  }

  last_updated_ = now + this_latency - original_latency_;
}

void DLPS::update_on_comm_end()
{
  XBT_DEBUG("Updating load of link '%s' on communication end", link_->get_cname());
  xbt_assert(is_enabled_,
             "Trying to update load of link '%s' while it is NOT enabled, aborting."
             " Please enable your link with sg_dlps_enable before trying to access any of its load metrics.",
             link_->get_cname());

  std::string link_name = link_->get_cname();
  double current_instantaneous_bytes_per_second = link_->get_usage();
  double now                                    = surf_get_clock();

  XBT_INFO("%s,finished,%.17f,%f\n", link_name.c_str(), now, current_instantaneous_bytes_per_second);

  update_load();
  last_busy_ = now;
}

s4u::Link* DLPS::get_s4u_link() {
  return link_;
}

bool DLPS::is_enabled() const
{
  return is_enabled_;
}

double DLPS::get_last_updated()
{
  return last_updated_;
}

double DLPS::get_last_busy()
{
  return last_busy_;
}

double DLPS::get_cumulated_bytes()
{
  update_load();
  return cumulated_bytes_;
}

double DLPS::get_min_bytes_per_second()
{
  update_load();
  return min_bytes_per_second_;
}

double DLPS::get_max_bytes_per_second()
{
  update_load();
  return max_bytes_per_second_;
}

double DLPS::get_average_bytes()
{
  update_load();

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
static void on_communicate(const simgrid::kernel::resource::NetworkAction& action)
{
  // XBT_DEBUG("on_communicate is called");
  for (auto* link : action.get_links()) {
    if (link == nullptr || link->get_sharing_policy() == simgrid::s4u::Link::SharingPolicy::WIFI)
      continue;

    auto dlps = link->get_iface()->extension<DLPS>();
    if (dlps->is_enabled()) {
      dlps->update_on_comm_start();
    }
  }
}

static void on_communication_state_change(const simgrid::kernel::resource::NetworkAction& action,
         simgrid::kernel::resource::Action::State) {
  for (auto const* link : action.get_links()) {
    if (link != nullptr && link->get_sharing_policy() != simgrid::s4u::Link::SharingPolicy::WIFI) {
      auto dlps = link->get_iface()->extension<DLPS>();
      if (dlps->is_enabled()) {
        if (action.get_state() == simgrid::kernel::resource::Action::State::FINISHED)
          dlps->update_on_comm_end();
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
