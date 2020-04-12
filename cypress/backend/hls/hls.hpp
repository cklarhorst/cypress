#ifndef CYPRESS_BACKEND_HLS_HPP
#define CYPRESS_BACKEND_HLS_HPP

#include <cypress/backend/hls/hls.hpp>
#include <cypress/core/backend.hpp>
#include <cypress/core/network_base_objects.hpp>
#include <cypress/util/json.hpp>

namespace cypress {

class HLSBackend : public Backend {
private:
	static std::string execHelper(const char* cmd);
	static std::string print_ident(int ident);
	static std::string change_ident(std::stringstream& inputStream, int ident);
	Logger *logger;
protected:
	void do_run(NetworkBase &network, Real duration) const override;
public:
	
	void create_source_population(const PopulationBase &pop, const int id, std::stringstream &declarations) const;
	void create_homogeneous_pop(const PopulationBase &pop, const int id,  std::stringstream &declarations, std::stringstream &neuron_updates_pre, std::stringstream &neuron_updates_post) const;
	void group_connect(const std::vector<PopulationBase> &populations, const ConnectionDescriptor &conn, const Real timestep, std::vector<std::vector<std::stringstream> > &spike_handling) const;

	/**
	 * Constructor of the PyNN backend. Throws an exception if the given PyNN
	 * backend does not exist.
	 *
	 * @param simulator is the name of the simulator backend to be used by PyNN.
	 * Use the static backends method to list available backends.
	 * @param setup contains additional setup information that should be passed
	 * to the backend.
	 */
	HLSBackend(const std::string &simulator, const Json &setup = Json());

	/**
	 * Destructor of the PyNN class.
	 */
	~HLSBackend() override;
	
	/**
	 * Returns the canonical name of the backend.
	 */
	std::string name() const override { return "test"; }
	
	/**
	 * Returns a set of neuron types which are supported by this backend. Trying
	 * to execute a network with other neurons than the ones specified in the
	 * result of this function will result in an exception.
	 *
	 * @return a set of neuron types supported by this particular backend
	 * instance.
	 */
	std::unordered_set<const NeuronType *> supported_neuron_types()
	    const override;
};

}

#endif /* CYPRESS_BACKEND_HLS_HPP */
