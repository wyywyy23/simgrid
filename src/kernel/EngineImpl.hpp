/* Copyright (c) 2016-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SIMGRID_KERNEL_ENGINEIMPL_HPP
#define SIMGRID_KERNEL_ENGINEIMPL_HPP

#include <simgrid/kernel/resource/Model.hpp>
#include <simgrid/s4u/Engine.hpp>
#include <simgrid/s4u/NetZone.hpp>
#include <simgrid/simix.hpp>

#include <map>
#include <string>
#include <unordered_map>

namespace simgrid {
namespace kernel {

class EngineImpl {
  std::map<std::string, s4u::Host*, std::less<>> hosts_;
  std::map<std::string, resource::LinkImpl*, std::less<>> links_;
  std::unordered_map<std::string, routing::NetPoint*> netpoints_;
  std::unordered_map<std::string, actor::ActorCodeFactory> registered_functions; // Maps function names to actor code
  actor::ActorCodeFactory default_function; // Function to use as a fallback when the provided name matches nothing
  std::vector<std::shared_ptr<resource::Model>> models_;
  std::unordered_map<resource::Model::Type, std::vector<resource::Model*>> models_by_type_;

  friend s4u::Engine;

public:
  EngineImpl() = default;

  EngineImpl(const EngineImpl&) = delete;
  EngineImpl& operator=(const EngineImpl&) = delete;
  virtual ~EngineImpl();

  void load_deployment(const std::string& file) const;
  void register_function(const std::string& name, const actor::ActorCodeFactory& code);
  void register_default(const actor::ActorCodeFactory& code);

  /**
   * @brief Add a model to engine list
   *
   * @param type Model type (network, disk, etc)
   * @param model Pointer to model
   * @param is_default Is this the default model for this type of resource in this exp
   */
  void add_model(resource::Model::Type type, std::shared_ptr<resource::Model> model, bool is_default = false);
  /** @brief Get current default model for a resource type */
  resource::Model* get_default_model(resource::Model::Type type) const;

  /** @brief Get list of models created for a resource type */
  const std::vector<resource::Model*>& get_model_list(resource::Model::Type type);

  /** @brief Get list of all models managed by this engine */
  const std::vector<std::shared_ptr<resource::Model>>& get_all_models() const { return models_; }

  routing::NetZoneImpl* netzone_root_ = nullptr;
  static EngineImpl* get_instance() { return simgrid::s4u::Engine::get_instance()->pimpl; }
  actor::ActorCodeFactory get_function(const std::string& name)
  {
    auto res = registered_functions.find(name);
    if (res == registered_functions.end())
      return default_function;
    else
      return res->second;
  }
};

} // namespace kernel
} // namespace simgrid

#endif
