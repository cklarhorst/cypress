#include <unistd.h>

#include <cypress/backend/resources.hpp>
#include <cypress/backend/hls/hls.hpp>
#include <cypress/core/network_base.hpp>
#include <cypress/core/network_base_objects.hpp>
#include <cypress/util/logger.hpp>
#include <sstream>
#include <iostream>
#include <fstream>

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

#define RECORD_EVERYTHING 1 //

namespace cypress {

/**
 * Helper to print vectors
 */
template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	out << "{";
	size_t last = v.size() - 1;
	for(size_t i = 0; i < v.size(); ++i) {
		out << (int)v[i];
		if (i != last) 
			out << ", ";
	}
	out << "}";
	return out;
}

/**
 * Idention helper
 */
std::string HLSBackend::print_ident(int ident)
{
	std::stringstream tmp;
	for (int i=0;i<ident;i++)
		tmp << "  ";
	return tmp.str();
}

/**
 * Idention helper
 */
std::string HLSBackend::change_ident(std::stringstream& inputStream, int ident) 
{ 
	std::string out;
	std::string temp; 
	
	out = "";
  
	while (getline(inputStream, temp, '\n')) { 
		out = out + print_ident(ident) + temp + '\n'; 
	} 
	return out; 
} 

/**
 * Helper to execute commands
 */
std::string HLSBackend::execHelper(const char* cmd)
{
	std::array<char, 128> buffer;
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
	if (!pipe) {
		throw std::runtime_error("popen() failed!");
	}

	while (1) 
	{
		size_t num_bytes = fread(buffer.data(),1,buffer.size(), pipe.get()); //todo: rewrite to parse input directly
		if (0==num_bytes)
			break;
		result.append(buffer.data(),num_bytes);
	}

	return result;
}

/**
 * Map containing the supported neuron types per simulator.
 */
static const std::unordered_map<std::string,
								std::unordered_set<const NeuronType *>>
	SUPPORTED_NEURON_TYPE_MAP = {
		{"__default__",
		{&SpikeSourceArray::inst(), &IfCondExp::inst(),
		&EifCondExpIsfaIsta::inst(), &IfCurrExp::inst()}}};

/**
 * Map containing connections
 */
static const std::unordered_map<std::string, std::string>
	SUPPORTED_CONNECTIONS = {
		{"AllToAllConnector", "AllToAllConnector"}
		};

HLSBackend::HLSBackend(__attribute__((unused)) const std::string &simulator, __attribute__((unused)) const Json &setup)
: logger(&global_logger())
{
}

HLSBackend::~HLSBackend() = default;

std::unordered_set<const NeuronType *> HLSBackend::supported_neuron_types() const
{
	return SUPPORTED_NEURON_TYPE_MAP.find("__default__")->second;
}

void HLSBackend::create_source_population(const PopulationBase &pop, const int id, std::stringstream &declarations) const
{
	declarations << "NEURON_SRC_T source_population_" << id << "[" << pop.size() << "][2000] = " << std::endl << "{" << std::endl; 
	size_t last = pop.size() - 1;
	for (size_t i = 0; i < pop.size(); ++i) {
		const auto neuron = pop[i];
		const std::vector<Real> &temp = neuron.parameters().parameters();
		declarations << temp;
		if (i != last)
			declarations << ", " << std::endl;
	}
	declarations << "};" << std::endl;
	declarations << "NEURON_SRC_ARRAY_POS_T source_population_" << id << "_pos[" << pop.size() << "] = {};" << std::endl;
}

void HLSBackend::create_homogeneous_pop(const PopulationBase &pop, const int id,  std::stringstream &declarations, std::stringstream &neuron_updates_pre, std::stringstream &neuron_updates_post) const
{
	const auto &params = pop[0].parameters();
	const auto &param_names = pop.type().parameter_names;

	declarations << "NEURON_CNT_T population_" << id << "_counter[" << pop.size() << "] = {};" << std::endl;

	int offset = abs(params[7]);

	for (size_t j = 0; j < param_names.size(); j++) {
		if (j == 5 || j == 6 || j == 7) 
			declarations << "NEURON_PARAM_T population_" << id << "_" << param_names[j] << " = " << params[j]+offset << ";" << std::endl;
		else
			declarations << "NEURON_PARAM_T population_" << id << "_" << param_names[j] << " = " << params[j] << ";" << std::endl;
	}
	for (size_t j = 0; j < pop.size(); j++) {
		neuron_updates_post << "if (population_" << id << "_counter[" << j << "] >= population_" << id << "_v_thresh) {" << std::endl;;
		auto neuron_inst = pop[j];
		if (neuron_inst.signals().is_recording(j) || RECORD_EVERYTHING) {
			neuron_updates_post << "#ifdef BIN_OUTPUT" << std::endl;
			neuron_updates_post << "    const unsigned char population_id = " << id << ";" << std::endl;
			neuron_updates_post << "    const unsigned char neuron_id     = " << j << ";" << std::endl;
			neuron_updates_post << "    std::cout.write((const char*)&population_id, 1);" << std::endl;
			neuron_updates_post << "    std::cout.write((const char*)&neuron_id, 1);" << std::endl;
			neuron_updates_post << "    std::cout.write((const char*)&time, 4);" << std::endl;
			neuron_updates_post << "#else" << std::endl;
			neuron_updates_post << "    std::cout << \"" << id << "," << j << ":\" << time << std::endl;" << std::endl;
			neuron_updates_post << "#endif" << std::endl;
		} else {
			neuron_updates_post << "// neuron " << j << " not marked for recording" << std::endl;
		}
		neuron_updates_post << "  population_" << id << "_counter[" << j << "] = population_" << id << "_v_rest;" << std::endl;
		neuron_updates_post << "}" << std::endl;
	}
	neuron_updates_pre  << "//Apply Leakage" << std::endl;
}


void HLSBackend::group_connect(__attribute__((unused)) const std::vector<PopulationBase> &populations,
								const ConnectionDescriptor &conn,
								__attribute__((unused)) const Real timestep, //TODO: impl?
								std::vector<std::vector<std::stringstream> > &spike_handling) const
{
	std::string conn_name = SUPPORTED_CONNECTIONS.find(conn.connector().name())->second;
	//const auto &params = conn.connector().synapse()->parameters();
	std::string name = conn.connector().synapse()->name();
	for (size_t s_id = conn.nid_src0(); s_id < (size_t) conn.nid_src1(); s_id++) {
		while (spike_handling[conn.pid_src()].size() <= s_id) {
			spike_handling[conn.pid_src()].push_back(std::stringstream());
		}
		for (size_t t_id = conn.nid_tar0(); t_id < (size_t) conn.nid_tar1(); t_id++) {
			std::string name = conn.connector().synapse()->name();
			//std::cout << params << std::endl;
			//params[0] // ex/in
			//params[1] // delay
			spike_handling[conn.pid_src()][s_id] << "population_" << conn.pid_tar() << "_counter[" << t_id << "] += population_" << conn.pid_tar() << "_tau_syn_E;" << std::endl;
		}
	}
}

void HLSBackend::do_run(NetworkBase &source, __attribute__((unused)) Real duration) const
{
	const std::vector<PopulationBase> &populations = source.populations();

	std::stringstream declarations;
	std::stringstream neuron_updates_pre;
	std::stringstream neuron_updates_post;

	for (size_t i = 0; i < populations.size(); i++) {
		if (populations.size() == 0) {
			logger->debug("HLS", "found empty");
			continue;
		}

		if (&populations[i].type() == &SpikeSourceArray::inst()) {
			logger->debug("HLS", "found source array");
			this->create_source_population(populations[i], i, declarations);
		}
		else {
			logger->debug("HLS", "found population");
			//bool homogeneous = populations[i].homogeneous_parameters(); todo: impl
			this->create_homogeneous_pop(populations[i], i, declarations, neuron_updates_pre, neuron_updates_post);
		}
	}
	std::vector< std::vector<std::stringstream> > spike_handling(populations.size());
	Real timestep = 0.1; //in ms;
	for (size_t i = 0; i < source.connections().size(); i++) {
		const auto &conn = source.connections()[i];
		auto it = SUPPORTED_CONNECTIONS.find(conn.connector().name());

		if (it != SUPPORTED_CONNECTIONS.end() &&
			conn.connector().group_connect(conn)) {
			// Group connections
			logger->debug("HLS", "found group connection");
			this->group_connect(populations, conn, timestep, spike_handling);
		}
		else {
			logger->error("HLS", "found other conntection -> not impl");

		}
	}
	// Code gen
	int ident = 0;
	std::ofstream output("hls_main.cpp", std::ios::out | std::ios::trunc);
	output << print_ident(ident) <<"#include <iostream>" << std::endl;
	output << print_ident(ident) << "#include <sstream>" << std::endl;
	output << print_ident(ident) << "#include <unistd.h>" << std::endl;
	output << print_ident(ident) << "// #define HLS" << std::endl;
	output << print_ident(ident) << "#define BIN_OUTPUT" << std::endl;
	output << print_ident(ident++) << "#ifdef HLS" << std::endl;
	output << print_ident(ident) << "#define NEURON_CNT_T ap_uint<6> " << std::endl;
	output << print_ident(ident) << "#define NEURON_SRC_T ap_uint<11> " << std::endl;
	output << print_ident(ident) << "#define NEURON_SRC_ARRAY_POS_T ap_uint<11> " << std::endl;
	output << print_ident(ident) << "#define NEURON_PARAM_T const ap_uint<5> " << std::endl;
	output << print_ident(ident-1) << "#else" << std::endl;
	output << print_ident(ident) << "#define NEURON_CNT_T uint32_t " << std::endl;
	output << print_ident(ident) << "#define NEURON_SRC_T uint32_t " << std::endl;
	output << print_ident(ident) << "#define NEURON_SRC_ARRAY_POS_T uint32_t " << std::endl;
	output << print_ident(ident) << "#define NEURON_PARAM_T const uint32_t " << std::endl;
	output << print_ident(--ident) << "#endif" << std::endl;
	output << print_ident(ident) << declarations.str();
	output << print_ident(ident++) << "int main()" << std::endl;
	output << print_ident(ident) << "{" << std::endl;
	output << print_ident(ident++) << "for (uint32_t time=0;time<1500;time++)" << std::endl; //TODO: Handle simulation time
	output << print_ident(ident) << "{" << std::endl;
	output << change_ident(neuron_updates_pre, ident) << std::endl;
	for (size_t p=0;p<populations.size();p++) {
		if (&populations[p].type() == &SpikeSourceArray::inst()) {
			for (size_t n=0;n<spike_handling[p].size();n++) {
				output << print_ident(ident++) << "if (source_population_" << p << "[" << n << "][source_population_" << p << "_pos" << "[" << n << "]] == time) {" << std::endl;
				output << change_ident(spike_handling[p][n], ident);
				output << print_ident(ident) << "source_population_" << p << "_pos[" << n << "]++;" << std::endl;
				output << print_ident(--ident) << "}" << std::endl;
			}
		}
	}
	output << change_ident(neuron_updates_post, ident) << std::endl;
	output << print_ident(--ident) << "}" << std::endl;
	output << print_ident(--ident) << "}" << std::endl;
	output.close();
	// Compiling
	std::string compiler_output = execHelper("g++ hls_main.cpp -Wall -o hls_main");
	logger->info("HLS", "Compiler output:");
	logger->info("HLS", compiler_output);
	logger->info("HLS", "Compiler output end");
	// Running
	std::string buf = execHelper("./hls_main");
	logger->debug("HLS", "Start parsing");
	// Parse results
	std::istringstream bufStream(buf);
	
	std::vector<std::vector<std::vector<uint32_t> > > pop_neur_time_buf;
	while (!bufStream.eof()) {
		uint8_t population_id;
		uint8_t neuron_id;
		uint32_t time;
		bufStream.read((char*)&population_id,1);
		bufStream.read((char*)&neuron_id,1);
		bufStream.read((char*)&time,4);
		
		if (pop_neur_time_buf.size()<=population_id)
			pop_neur_time_buf.resize(population_id+1);
		if (pop_neur_time_buf[population_id].size()<=neuron_id)
			pop_neur_time_buf[population_id].resize(neuron_id+1);
		pop_neur_time_buf[population_id][neuron_id].push_back(time);
	}
	
	logger->info("HLS", "Parsed " + std::to_string(pop_neur_time_buf.size()) + " #populations");
	
	for (size_t population_id=0;population_id<pop_neur_time_buf.size();population_id++)
	{
		for (size_t neuron_id=0;neuron_id<pop_neur_time_buf[population_id].size();neuron_id++)
		{
			logger->info("HLS", "Parsed " + std::to_string(pop_neur_time_buf[population_id].size()) + " #neurons");
			logger->info("HLS", "Parsed " + std::to_string(pop_neur_time_buf[population_id][neuron_id].size()) + " #times");
			if (pop_neur_time_buf[population_id][neuron_id].size()==0)
				continue;
			auto data = std::make_shared<Matrix<Real>>(pop_neur_time_buf[population_id][neuron_id].size(), 1);
			for (size_t j = 0; j < pop_neur_time_buf[population_id][neuron_id].size(); j++) {
				(*data)(j, 0) = pop_neur_time_buf[population_id][neuron_id][j];
			}
			auto &pop = populations[population_id];
			auto neuron_inst = pop[neuron_id];
			neuron_inst.signals().data(0, std::move(data));
		}
	}
	logger->debug("HLS", "End do_run");
}

}
