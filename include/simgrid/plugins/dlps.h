/* Copyright (c) 2017-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SIMGRID_PLUGINS_DLPS_H_
#define SIMGRID_PLUGINS_DLPS_H_

#include <simgrid/config.h>
#include <simgrid/forward.h>
#include <xbt/base.h>

SG_BEGIN_DECL

XBT_PUBLIC void sg_dlps_plugin_init();
XBT_PUBLIC void sg_dlps_enable(const_sg_link_t link);
XBT_PUBLIC void sg_dlps_enable(const_sg_link_t link);
XBT_PUBLIC void sg_dlps_reset(const_sg_link_t link);
XBT_PUBLIC double sg_dlps_get_cum_load(const_sg_link_t link);
XBT_PUBLIC double sg_dlps_get_avg_load(const_sg_link_t link);
XBT_PUBLIC double sg_dlps_get_min_instantaneous_load(const_sg_link_t link);
XBT_PUBLIC double sg_dlps_get_max_instantaneous_load(const_sg_link_t link);

SG_END_DECL

#endif
