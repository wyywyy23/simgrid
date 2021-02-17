/* Copyright g(c) 2019-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/s4u.hpp"
#include "xbt/asserts.h"
#include "xbt/log.h"

XBT_LOG_NEW_DEFAULT_CATEGORY(ptask_L07_tes, "[usage] ptask_LO7 <platform-file>");

namespace sg4 = simgrid::s4u;

/* We need a separate actor so that it can sleep after each test */
static void main_dispatcher()
{
  sg4::Engine* e = sg4::Engine::get_instance();
  double start_time;
  double end_time;
  std::vector<sg4::Host*> hosts = e->get_all_hosts();

  XBT_INFO("TEST: Create and run a sequential execution.");
  XBT_INFO("------------------------------------------------------------");
  XBT_INFO("Have to compute 1 flop on a 1 flop/s host.");
  XBT_INFO("Should be done in exactly one second.");
  start_time = e->get_clock();
  sg4::Exec::init()->set_flops_amount(1)->set_host(hosts[0])->wait();
  end_time = e->get_clock();
  XBT_INFO("Actual result: computing 1 flop at 1 flop/s takes %.2f seconds.", end_time - start_time);
  XBT_INFO("\n");

  sg4::this_actor::sleep_for(5);

  XBT_INFO("TEST: Create and run a parallel execution on 2 homogeneous hosts.");
  XBT_INFO("------------------------------------------------------------");
  XBT_INFO("Have to compute 2 flops across two hosts running at 1 flop/s.");
  XBT_INFO("Should be done in exactly one second.");
  start_time = e->get_clock();
  sg4::Exec::init()->set_flops_amounts(std::vector<double>({1.0, 1.0}))
                   ->set_hosts(std::vector<sg4::Host*>({hosts[0], hosts[1]}))
                   ->wait();
  end_time = e->get_clock();
  XBT_INFO("Actual result: computing 2 flops on 2 hosts at 1 flop/s takes %.2f seconds.", end_time - start_time);
  XBT_INFO("\n");

  sg4::this_actor::sleep_for(5);

  XBT_INFO("TEST: Create and run a parallel execution on 2 heterogeneous hosts.");
  XBT_INFO("------------------------------------------------------------");
  XBT_INFO("Have to compute 2 flops across two hosts, one running at 1 flop/s and one at 2 flop/s.");
  XBT_INFO("Should be done in exactly one second.");
  start_time = e->get_clock();
  sg4::Exec::init()->set_flops_amounts(std::vector<double>({1.0, 1.0}))
                   ->set_hosts(std::vector<sg4::Host*>({hosts[1], hosts[2]}))
                   ->wait();
  end_time = e->get_clock();
  XBT_INFO("Actual result: computing 2 flops on 2 heterogeneous hosts takes %.2f seconds.", end_time - start_time);
  XBT_INFO("\n");

  sg4::this_actor::sleep_for(5);

  XBT_INFO("TEST: Latency test between hosts connected by a shared link.");
  XBT_INFO("------------------------------------------------------------");
  XBT_INFO("Have to send 1B from one host to another at 1Bps with a latency of 500ms.");
  XBT_INFO("Should be done in 1.5 seconds (500ms latency + 1s transfert).");
  start_time = e->get_clock();
  sg4::Comm::sendto_async(hosts[0], hosts[1], 1.0)->wait();
  end_time = e->get_clock();
  XBT_INFO("Actual result: sending 1 byte on a shared link at 1Bps + 500ms takes %.2f seconds.", end_time - start_time);
  XBT_INFO("\n");

  sg4::this_actor::sleep_for(5);

  XBT_INFO("TEST: Latency test between hosts connected by a fatpipe link.");
  XBT_INFO("------------------------------------------------------------");
  XBT_INFO("Have to send 1B from one host to another at 1Bps with a latency of 500ms.");
  XBT_INFO("Should be done in 1.5 seconds (500ms latency + 1s transfert).");
  start_time = e->get_clock();
  sg4::Comm::sendto_async(hosts[0], hosts[2], 1.0)->wait();
  end_time = e->get_clock();
  XBT_INFO("Actual result: sending 1 byte on a fatpipe link at 1Bps + 500ms takes %.2f seconds.", end_time - start_time);
  XBT_INFO("\n");

  sg4::this_actor::sleep_for(5);

  XBT_INFO("TEST: Latency test between hosts connected by a 3-link route.");
  XBT_INFO("------------------------------------------------------------");
  XBT_INFO("Have to send 1B from one host to another at 1Bps with a latency of 2 x 500ms + 1s.");
  XBT_INFO("Should be done in 3 seconds (2 x 500ms + 1s latency + 1s transfert).");
  start_time = e->get_clock();
  sg4::Comm::sendto_async(hosts[1], hosts[2], 1.0)->wait();
  end_time = e->get_clock();
  XBT_INFO("Actual result: sending 1 byte on a 3-link route at 1Bps + 2,500ms takes %.2f seconds.", end_time - start_time);
  XBT_INFO("\n");

  sg4::this_actor::sleep_for(5);

   XBT_INFO("TEST: Latency test between hosts connected by a shared link with 2 comms in same direction.");
   XBT_INFO("------------------------------------------------------------");
   XBT_INFO("Have to send 2 x 1B from one host to another at 1Bps with a latency of 500ms.");
   XBT_INFO("Should be done in 2.5 seconds (500ms latency + 2s transfert).");
   start_time = e->get_clock();
   sg4::CommPtr c1 = sg4::Comm::sendto_async(hosts[0], hosts[1], 1.0);
   sg4::CommPtr c2 = sg4::Comm::sendto_async(hosts[0], hosts[1], 1.0);
   c1->wait();
   c2->wait();
   end_time = e->get_clock();
   XBT_INFO("Actual result: sending 2x1 bytes on a shared link at 1Bps + 500ms takes %.2f seconds.", end_time - start_time);
   XBT_INFO("\n");

   sg4::this_actor::sleep_for(5);

   XBT_INFO("TEST: Latency test between hosts connected by a fatpipe link with 2 comms in same direction.");
   XBT_INFO("------------------------------------------------------------");
   XBT_INFO("Have to send 2 x 1B from one host to another at 1Bps with a latency of 500ms.");
   XBT_INFO("Should be done in 1.5 seconds (500ms latency + 1s transfert).");
   start_time = e->get_clock();
   c1 = sg4::Comm::sendto_async(hosts[0], hosts[2], 1.0);
   c2 = sg4::Comm::sendto_async(hosts[0], hosts[2], 1.0);
   c1->wait();
   c2->wait();
   end_time = e->get_clock();
   XBT_INFO("Actual result: sending 2x1 bytes on a fatpipe link at 1Bps + 500ms takes %.2f seconds.", end_time - start_time);
   XBT_INFO("\n");

   sg4::this_actor::sleep_for(5);

   XBT_INFO("TEST: Latency test between hosts connected by a 3-link route with 2 comms in same direction.");
   XBT_INFO("------------------------------------------------------------");
   XBT_INFO("Have to send 2 x 1B from one host to another at 1Bps with a latency of 2 x 500ms + 1s.");
   XBT_INFO("Should be done in 4 seconds (2 x 500ms + 1s latency + 2s transfert).");
   start_time = e->get_clock();
   c1 = sg4::Comm::sendto_async(hosts[1], hosts[2], 1.0);
   c2 = sg4::Comm::sendto_async(hosts[1], hosts[2], 1.0);
   c1->wait();
   c2->wait();
   end_time = e->get_clock();
   XBT_INFO("Actual result: sending 2x1 bytes on a 3-link route at 1Bps + 2,500ms takes %.2f seconds.",
            end_time - start_time);
   XBT_INFO("\n");

   sg4::this_actor::sleep_for(5);

   XBT_INFO("TEST: Latency test between hosts connected by a shared link with 2 comms in opposite direction.");
   XBT_INFO("------------------------------------------------------------");
   XBT_INFO("Have to send 1B between two hosts in each direction at 1Bps with a latency of 500ms.");
   XBT_INFO("Should be done in 2.5 seconds (500ms latency + 2s transfert).");
   start_time = e->get_clock();
   c1 = sg4::Comm::sendto_async(hosts[0], hosts[1], 1.0);
   c2 = sg4::Comm::sendto_async(hosts[1], hosts[0], 1.0);
   c1->wait();
   c2->wait();
   end_time = e->get_clock();
   XBT_INFO("Actual result: sending 1 byte in both directions on a shared link at 1Bps + 500ms takes %.2f seconds.",
            end_time - start_time);
   XBT_INFO("\n");

   sg4::this_actor::sleep_for(5);

   XBT_INFO("TEST: Latency test between hosts connected by a fatpipe link with 2 comms in opposite direction.");
   XBT_INFO("------------------------------------------------------------");
   XBT_INFO("Have to send 1B between two hosts in each direction at 1Bps with a latency of 500ms.");
   XBT_INFO("Should be done in 1.5 seconds (500ms latency + 1s transfert).");
   start_time = e->get_clock();
   c1 = sg4::Comm::sendto_async(hosts[0], hosts[2], 1.0);
   c2 = sg4::Comm::sendto_async(hosts[2], hosts[0], 1.0);
   c1->wait();
   c2->wait();
   end_time = e->get_clock();
   XBT_INFO("Actual result: sending 1 byte in both directions on a fatpipe link at 1Bps + 500ms takes %.2f seconds.",
            end_time - start_time);
   XBT_INFO("\n");

   sg4::this_actor::sleep_for(5);

   XBT_INFO("TEST: Latency test between hosts connected by a 3-link route with 2 comms in opposite direction.");
   XBT_INFO("------------------------------------------------------------");
   XBT_INFO("Have to send 1B between two hosts in each direction at 1Bps with a latency of 2 x 500ms + 1s.");
   XBT_INFO("Should be done in 4 seconds (2 x 500ms + 1s latency + 2s transfert).");
   start_time = e->get_clock();
   c1 = sg4::Comm::sendto_async(hosts[1], hosts[2], 1.0);
   c2 = sg4::Comm::sendto_async(hosts[2], hosts[1], 1.0);
   c1->wait();
   c2->wait();
   end_time = e->get_clock();
   XBT_INFO("Actual result: sending 1 byte in both directions on a 3-link route at 1Bps + 2,500ms takes %.2f seconds.",
             end_time - start_time);
   XBT_INFO("\n");
}

int main(int argc, char** argv)
{
  simgrid::s4u::Engine engine(&argc, argv);
  engine.load_platform(argv[1]);
  simgrid::s4u::Actor::create("dispatcher", simgrid::s4u::Host::by_name("cpu0"), main_dispatcher);
  engine.run();

  return 0;
}

