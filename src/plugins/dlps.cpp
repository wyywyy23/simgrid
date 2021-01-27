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
  void update();
  double get_average_bytes();

  /// Getter methods.
  bool is_enabled() const;
  s4u::Link* get_s4u_link();
  double get_cumulated_bytes();
  double get_min_bytes_per_second();
  double get_max_bytes_per_second();

private:
  s4u::Link* link_{};      /*< The link onto which this data is enabled*/
  bool is_enabled_{false}; /*<Whether the link is enabled or not*/

  double cumulated_bytes_{};      /*< Cumulated load since last reset*/
  double min_bytes_per_second_{}; /*< Minimum instantaneous load observed since last reset*/
  double max_bytes_per_second_{}; /*< Maximum instantaneous load observed since last reset*/
  double last_reset_{};           /*< Timestamp of the last reset (init timestamp by default)*/
  double last_updated_{};         /*< Timestamp of the last energy update event*/
};

xbt::Extension<s4u::Link, DLPS> DLPS::EXTENSION_ID;

DLPS::DLPS(simgrid::s4u::Link* ptr) : link_(ptr), is_enabled_(false)
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

void DLPS::update()
{
  XBT_DEBUG("Updating load of link '%s'", link_->get_cname());
  xbt_assert(is_enabled_,
             "Trying to update load of link '%s' while it is NOT enabled, aborting."
             " Please enable your link with sg_dlps_enable before trying to access any of its load metrics.",
             link_->get_cname());

  std::string link_name = link_->get_cname();
  double current_instantaneous_bytes_per_second = link_->get_usage();
  double now                                    = surf_get_clock();
  
//  if (link_name[1] != 'a' and link_name[1] != 'o') {
//    XBT_INFO("%s,usage,%.17f,%f\n", link_->get_cname(), now, current_instantaneous_bytes_per_second);
//  }

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

s4u::Link* DLPS::get_s4u_link() {
  return link_;
}

bool DLPS::is_enabled() const
{
  return is_enabled_;
}
double DLPS::get_cumulated_bytes()
{
  update();
  return cumulated_bytes_;
}
double DLPS::get_min_bytes_per_second()
{
  update();
  return min_bytes_per_second_;
}
double DLPS::get_max_bytes_per_second()
{
  update();
  return max_bytes_per_second_;
}

double DLPS::get_average_bytes()
{
  update();

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
      dlps->update();

      std::string link_name = dlps->get_s4u_link()->get_cname();
      double current_instantaneous_bytes_per_second = dlps->get_s4u_link()->get_usage();
      double now                                    = surf_get_clock();

      if (link_name[1] != 'a' and link_name[1] != 'o') {
        XBT_INFO("%s,on_communicate,%.17f,%f\n", link_name.c_str(), now, current_instantaneous_bytes_per_second);
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
  simgrid::s4u::Link::on_communicate.connect(&on_communicate);
  simgrid::s4u::Link::on_state_change.connect([](simgrid::s4u::Link const& link) {
    if (link.get_sharing_policy() != simgrid::s4u::Link::SharingPolicy::WIFI) {
      auto dlps = link.extension<DLPS>();
      if (dlps->is_enabled())
        dlps->update();
    }
  });
  simgrid::s4u::Link::on_communication_state_change.connect(
      [](simgrid::kernel::resource::NetworkAction const& action,
         simgrid::kernel::resource::Action::State) {
        for (auto const* link : action.get_links()) {
          if (link != nullptr && link->get_sharing_policy() != simgrid::s4u::Link::SharingPolicy::WIFI) {
            auto dlps = link->get_iface()->extension<DLPS>();
            if (dlps->is_enabled()) {
              dlps->update();

              std::string link_name = dlps->get_s4u_link()->get_cname();
              double current_instantaneous_bytes_per_second = dlps->get_s4u_link()->get_usage();
              double now                                    = surf_get_clock();
	      std::string current_state                     = "";

	      switch (action.get_state()) {
		case simgrid::kernel::resource::Action::State::INITED : 
			current_state = "initialized"; break;
		case simgrid::kernel::resource::Action::State::STARTED : 
			current_state = "started"; break;
		case simgrid::kernel::resource::Action::State::FAILED : 
			current_state = "failed"; break;
		case simgrid::kernel::resource::Action::State::FINISHED : 
			current_state = "finished"; break;
		case simgrid::kernel::resource::Action::State::IGNORED : 
			current_state = "ignored"; break;
	      }

              if (link_name[1] != 'a' and link_name[1] != 'o') {
                XBT_INFO("%s,%s,%.17f,%f\n", link_name.c_str(), current_state.c_str(), now, current_instantaneous_bytes_per_second);
              }
	    }
	      
          }
        }
      });
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
