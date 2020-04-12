/*
 *  Cypress -- C++ Spiking Neural Network Simulation Framework
 *  Copyright (C) 2016  Andreas Stöckel
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Include first to avoid "_POSIX_C_SOURCE redefined" warning
#include <cypress/backend/pynn/pynn.hpp>
#include <cypress/core/network_base.hpp>

#include <algorithm>
#include <memory>
#include <sstream>

#include <cypress/core/backend.hpp>
#include <cypress/core/data.hpp>
#include <cypress/core/exceptions.hpp>
#include <cypress/core/network_base_objects.hpp>
#include <cypress/core/transformation.hpp>

#include <cypress/backend/brainscales/brainscales_lib.hpp>
#include <cypress/backend/brainscales/slurm.hpp>
#include <cypress/backend/hls/hls.hpp>
#include <cypress/backend/nest/nest.hpp>
#include <cypress/backend/genn/genn_lib.hpp>
#include <cypress/backend/nmpi/nmpi.hpp>
#include <cypress/backend/serialize/to_json.hpp>

#include <cypress/transformations/registry.hpp>

#include <cypress/util/json.hpp>
#include <cypress/util/logger.hpp>

#include <dlfcn.h>

namespace cypress {
namespace internal {
/**
 * Spiking neural network data container. Stores all the data describing the
 * network.
 */
class NetworkData {
private:
	friend NetworkBase;

	/**
	 * Contains information about the runtime.
	 */
	NetworkRuntime m_runtime;

	/**
	 * Vector containing the PopulationData instances.
	 */
	std::vector<std::shared_ptr<PopulationData>> m_populations;

	/**
	 * Vector containing all connections.
	 */
	mutable std::vector<ConnectionDescriptor> m_connections;

	/**
	 * Flag indicating whether the connections are currently sorted.
	 */
	mutable bool m_connections_sorted;

public:
	/**
	 * Logger that is being used.
	 */
	Logger *logger;

	/**
	 * Flag indicating whether lossy transformations should be used when
	 * executing the network.
	 */
	bool use_lossy_trafos;

	/**
	 * Set of disabled transformation identifiers.
	 */
	std::unordered_set<std::string> disabled_trafo_ids;

	/**
	 * Default constructor of the NetworkData class.
	 */
	NetworkData()
	    : m_runtime({}),
	      m_connections_sorted(true),
	      logger(&global_logger()),
	      use_lossy_trafos(true){};

	/**
	 * Creates an independent NetworkData instance.
	 */
	NetworkData clone()
	{
		NetworkData res = *this;
		for (auto &sp : res.m_populations) {
			sp = std::make_shared<PopulationData>(*sp);
		}
		return res;
	}

	/**
	 * Returns a reference at the vector containing the population data.
	 */
	std::vector<std::shared_ptr<PopulationData>> &populations()
	{
		return m_populations;
	}

	/**
	 * Returns the indices for the populations which match the given search
	 * criteria.
	 */
	std::vector<PopulationIndex> populations(const std::string &name,
	                                         const NeuronType &type)
	{
		std::vector<PopulationIndex> res;
		for (size_t pid = 0; pid < m_populations.size(); pid++) {
			const PopulationData &pop = *m_populations[pid];
			if ((name.empty() || pop.name() == name) &&
			    (&type == &NullNeuron::inst() || pop.type() == &type)) {
				res.push_back(pid);
			}
		}
		return res;
	}

	/**
	 * Adds the given connector to the connection list.
	 */
	void connect(const ConnectionDescriptor &descr)
	{
		// Make sure the target population is not a spike source
		if (m_populations[descr.pid_tar()]->type()->spike_source) {
			throw InvalidConnectionException(
			    "Spike sources are not valid connection targets.");
		}

		// Assemble the connection descriptor and check its validity
		if (!descr.valid()) {
			throw InvalidConnectionException(
			    "The source and target population sizes do not match the size "
			    "expected by the chosen connector.");
		}

		// Append the descriptor to the connection list, update the sorted flag
		m_connections.emplace_back(std::move(descr));
		m_connections_sorted = (m_connections.size() <= 1) ||
		                       (m_connections_sorted &&
		                        m_connections[m_connections.size() - 2] <
		                            m_connections[m_connections.size() - 1]);
	}

	/**
	 * Returns the list of connections. Makes sure the returned connections are
	 * sorted.
	 */
	std::vector<ConnectionDescriptor> &connections()
	{
		if (!m_connections_sorted) {
			std::sort(m_connections.begin(), m_connections.end());
			m_connections_sorted = true;
		}
		return m_connections;
	}

	/**
	 * Returns the first connection with name
	 */
	ConnectionDescriptor &connections(std::string name)
	{
		for (size_t i = 0; i < m_connections.size(); i++) {
			if (m_connections[i].label() == name) {
				return m_connections[i];
			}
		}
		throw NoSuchPopulationException(std::string("Connection with name \"") +
		                                name + "\" does not exist");
	}
};
}  // namespace internal

/*
 * Class NetworkBase
 */

NetworkBase::NetworkBase() : m_impl(std::make_shared<internal::NetworkData>())
{
	// Do nothing here
}

NetworkBase::~NetworkBase() = default;

Logger &NetworkBase::logger() const { return *m_impl->logger; }

void NetworkBase::logger(Logger &logger) { m_impl->logger = &logger; }

NetworkBase NetworkBase::clone() const
{
	return NetworkBase(
	    std::make_shared<internal::NetworkData>(m_impl->clone()));
}

void NetworkBase::connect(PopulationIndex pid_src, NeuronIndex nid_src0,
                          NeuronIndex nid_src1, PopulationIndex pid_tar,
                          NeuronIndex nid_tar0, NeuronIndex nid_tar1,
                          std::unique_ptr<Connector> connector,
                          const char *name)
{
	m_impl->connect(ConnectionDescriptor(pid_src, nid_src0, nid_src1, pid_tar,
	                                     nid_tar0, nid_tar1,
	                                     std::move(connector), name));
}

PopulationIndex NetworkBase::create_population_index(
    size_t size, const NeuronType &type, const NeuronParameters &params,
    const NeuronSignals &signals, const std::string &name)
{
	// Create a new population data store
	auto data = std::make_shared<PopulationData>(size, &type, name);

	// Copy the given parameters to the new population
	NeuronParameters(data, 0, size) = params;
	NeuronSignals(data, 0, size) = signals;

	// Append the population the the existing list of populations
	m_impl->populations().emplace_back(data);
	return population_count() - 1;
}

size_t NetworkBase::population_count() const
{
	return m_impl->populations().size();
}

size_t NetworkBase::population_count(const NeuronType &type) const
{
	size_t res = 0;
	for (const auto &pop : m_impl->populations()) {
		if (pop->type() == &type) {
			res++;
		}
	}
	return res;
}

size_t NetworkBase::neuron_count() const
{
	size_t res = 0;
	for (const auto &pop : m_impl->populations()) {
		res += pop->size();
	}
	return res;
}

size_t NetworkBase::neuron_count(const NeuronType &type) const
{
	size_t res = 0;
	for (const auto &pop : m_impl->populations()) {
		if (pop->type() == &type) {
			res += pop->size();
		}
	}
	return res;
}

std::shared_ptr<PopulationData> NetworkBase::population_data(
    PopulationIndex pid)
{
	return m_impl->populations()[pid];
}

std::shared_ptr<const PopulationData> NetworkBase::population_data(
    PopulationIndex pid) const
{
	return m_impl->populations()[pid];
}

PopulationBase NetworkBase::population(const std::string &name)
{
	auto pops = populations(name);
	if (pops.empty()) {
		throw NoSuchPopulationException(std::string("Population with name \"") +
		                                name + "\" does not exist");
	}
	return pops.back();
}

PopulationBase NetworkBase::population(PopulationIndex pid)
{
	return PopulationBase(m_impl, pid);
}

const PopulationBase NetworkBase::population(PopulationIndex pid) const
{
	return PopulationBase(m_impl, pid);
}

PopulationBase NetworkBase::operator[](PopulationIndex pid)
{
	return population(pid);
}

const PopulationBase NetworkBase::operator[](PopulationIndex pid) const
{
	return population(pid);
}

const std::vector<PopulationBase> NetworkBase::populations(
    const std::string &name, const NeuronType &type) const
{
	std::vector<PopulationIndex> pids = m_impl->populations(name, type);
	std::vector<PopulationBase> res;
	for (PopulationIndex pid : pids) {
		res.emplace_back(*this, pid);
	}
	return res;
}

const std::vector<PopulationBase> NetworkBase::populations(
    const NeuronType &type) const
{
	return populations(std::string(), type);
}

const std::vector<ConnectionDescriptor> &NetworkBase::connections() const
{
	return m_impl->connections();
}

const ConnectionDescriptor &NetworkBase::connection(std::string name) const
{
	return m_impl->connections(name);
}

static std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	std::stringstream ss(s);
	std::string elem;
	while (getline(ss, elem, delim)) {
		elems.push_back(elem);
	}
	return elems;
}

static std::string join(const std::vector<std::string> &elems, char delim)
{
	bool first = true;
	std::stringstream ss;
	for (const std::string &elem : elems) {
		if (!first) {
			ss << delim;
		}
		first = false;
		ss << elem;
	}
	return ss.str();
}

void NetworkBase::update_connection(std::unique_ptr<Connector> connector,
                                    const char *name)
{
	std::vector<ConnectionDescriptor> &connections = m_impl->connections();
	int index = -1;
	for (size_t i = 0; i < connections.size(); i++) {
		if (connections[i].label() == name) {
			if (index > -1) {
				throw std::invalid_argument(
				    "The name of the connection is ambiguous");
			}
			index = i;
		}
	}
	if (index < 0) {
		throw std::invalid_argument("The name of the connection doesn't exist");
	}
	m_impl->connections()[index].update_connector(std::move(connector));
}

std::unique_ptr<Backend> NetworkBase::make_backend(std::string backend_id,
                                                   int argc, const char *argv[],
                                                   Json setup)
{
	// Split the input string at the "=" -- the lefthand part contains the
	// backend identifier, the right-hand part possible Json simulator setup
	// data.
	size_t config_pos = backend_id.find("=");
	if (config_pos != std::string::npos) {
		if (!setup.is_null()) {
			throw std::invalid_argument(
			    "Setup data present in the backend identifier, but explicit "
			    "setup given in the call to make_backend()!");
		}
		setup = Json::parse(backend_id.substr(
		    config_pos + 1, backend_id.size() - config_pos + 1));
		backend_id = backend_id.substr(0, config_pos);
	}
	std::vector<std::string> elems = split(backend_id, '.');

	// Make sure the backend identifier is not empty
	if (elems.empty()) {
		throw std::invalid_argument("Backend ID must not be empty!");
	}

	// Forward calls starting with "nmpi" to the NMPI platform
	if (elems[0] == "nmpi") {
		if (argc <= 0 || !argv) {
			throw std::invalid_argument(
			    "Command line arguments were not passed to run(), but these "
			    "are required for the NMPI backend!");
		}
		elems.erase(elems.begin());  // Remove the first element
		if (elems.empty()) {
			throw std::invalid_argument(
			    "Expected another backend name following \"nmpi\"!");
		}
		if (elems[0] == "nmpm1" || elems[0] == "BrainScaleS" ||
		    elems[0] == "brainscales") {

			if (!setup.is_null()) {
				elems[0] = "nmpm1=" + setup.dump();
			}
			else {
				elems[0] = "nmpm1";
			}
			return std::make_unique<NMPI>(join(elems, '.') + "=" + setup.dump(),
			                              argc, argv);
		}

		auto backend = make_backend(join(elems, '.'), argc, argv, setup);
		if (dynamic_cast<PyNN *>(backend.get()) == nullptr) {
			throw std::invalid_argument(
			    "NMPI backend only works in conjunction with PyNN "
			    "or BrainScaleS backends!");
		}
		std::unique_ptr<PyNN> pynn_backend(dynamic_cast<PyNN *>(backend.get()));
		backend.release();
		return std::make_unique<NMPI>(std::move(pynn_backend), argc, argv);
	}
	else if (elems[0] == "pynn") {
		elems.erase(elems.begin());  // Remove the first element
		if (elems.empty()) {
			throw std::invalid_argument(
			    "Expected another backend name following \"pynn\"!");
		}
		return std::make_unique<PyNN>(join(elems, '.'), setup);
	}
	else if (elems[0] == "slurm") {
		elems.erase(elems.begin());  // Remove the first element
		if (elems.empty()) {
			throw std::invalid_argument(
			    "Expected another backend name following \"slurm\"!");
		}
		return std::make_unique<Slurm>(join(elems, '.'), setup);
	}
	else if (elems[0] == "nmpm1") {
		return std::unique_ptr<Backend>(
		    BS_Lib::instance().create_bs_backend(setup));
	}
	else if (elems[0] == "nest") {
		return std::make_unique<NEST>(setup);
	}
	else if (elems[0] == "genn") {
		return std::unique_ptr<Backend>(
		    GENN_Lib::instance().create_genn_backend(setup));
	}
	else if (elems[0] == "json") {
		elems.erase(elems.begin());
		return std::make_unique<ToJson>(join(elems, '.'), setup);
	} 
	else if (elems[0] == "hls") {
		elems.erase(elems.begin());
		return std::make_unique<HLSBackend>(join(elems, '.'), setup);
	}
	else {
		return std::make_unique<PyNN>(join(elems, '.'), setup);
	}
	return nullptr;
}

bool NetworkBase::use_lossy_trafos() const { return m_impl->use_lossy_trafos; }

void NetworkBase::use_lossy_trafos(bool use_lossy) const
{
	m_impl->use_lossy_trafos = use_lossy;
}

const std::unordered_set<std::string> &NetworkBase::disabled_trafo_ids() const
{
	return m_impl->disabled_trafo_ids;
}

std::unordered_set<std::string> &NetworkBase::disabled_trafo_ids()
{
	return m_impl->disabled_trafo_ids;
}

void NetworkBase::run(const Backend &backend, Real duration)
{
	// Automatically deduce the duration if none was given
	if (duration <= 0) {
		duration = std::round(this->duration() + 1000.0);
	}

	// Run the network through the transformation machinery, make sure all
	// transformations are registered
	logger().info("cypress", "Executing network...");
	transformations::register_();
	Transformations::run(backend, *this, TransformationAuxData{duration},
	                     m_impl->disabled_trafo_ids, m_impl->use_lossy_trafos);

	// Print some execution summary
	auto rt = runtime();
	std::stringstream ss;
	ss << "Done. Execution took " << rt.total << "s (simulation " << rt.sim
	   << "s, initialization " << rt.initialize << "s, finalization "
	   << rt.finalize << "s)";
	logger().info("cypress", ss.str());
}

void NetworkBase::run(const std::string &backend_id, Real duration, int argc,
                      const char *argv[])
{
	auto backend = make_backend(backend_id, argc, argv);
	run(*backend, duration);
}

Real NetworkBase::duration() const
{
	Real res = 0.0;
	for (const auto &population : populations()) {
		if (&population.type() == &SpikeSourceArray::inst()) {
			const NeuronIndex nid_end =
			    population.homogeneous_parameters() ? 1 : population.size();
			for (NeuronIndex nid = 0; nid < nid_end; nid++) {
				auto &params = population[nid].parameters().parameters();
				if (params.size() > 0) {
					// Note: the spike times are supposed to be sorted!
					res = std::max(res, params[params.size() - 1]);
				}
			}
		}
	}
	return res;
}

NetworkRuntime NetworkBase::runtime() const { return m_impl->m_runtime; }

void NetworkBase::runtime(const NetworkRuntime &runtime)
{
	m_impl->m_runtime = runtime;
}
}  // namespace cypress
