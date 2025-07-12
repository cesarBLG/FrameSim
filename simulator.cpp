#include "simulator.h"
#include "nonsparsesim.h"
#include <iostream>
//#define RANDOMIZE_FLIPS
//#define CHECK_FT
// Flips an error for a given shot
// If there was an error, erase it
// If there was no error, insert it
void FrameSimulator::flip_error(size_t shot, int qubit, int type)
{
    auto it = errors.find(shot);
    if (it == errors.end())
        it = errors.emplace(shot, std::pair<std::set<int>,std::set<int>>()).first;
    flip_error(it, qubit, type);
}
// Flips an error for a given shot, using an iterator that points to the shot
// Advances iterator to next shot with an error
void FrameSimulator::flip_error(std::map<size_t,std::pair<std::set<int>,std::set<int>>>::iterator &it, int qubit, int type)
{
    if ((type & ERROR_X) && !it->second.first.erase(qubit))
        it->second.first.insert(qubit);
    if ((type & ERROR_Z) && !it->second.second.erase(qubit))
        it->second.second.insert(qubit);
    if (it->second.first.empty() && it->second.second.empty())
        it = errors.erase(it);
    else
        ++it;
}
// Erases an error for a given shot
// Advances iterator to next shot with error
void FrameSimulator::reset_error(std::map<size_t,std::pair<std::set<int>,std::set<int>>>::iterator &it, int qubit, int type)
{
    // FIXME!!! Not implemented for Y errors
    if (type == ERROR_X)
        it->second.first.erase(qubit);
    else if (type == ERROR_Z)
        it->second.second.erase(qubit);
    else
        abort();
    if (it->second.first.empty() && it->second.second.empty())
        it = errors.erase(it);
    else
        ++it;
}
// Hadamard gate
// Exchanges X and Z errors in the target qubit
void FrameSimulator::h(int qubit)
{
    for (auto &[shot, err] : errors) {
        auto it1 = err.first.find(qubit);
        auto it2 = err.second.find(qubit);
        if (it1 != err.first.end()) {
            if (it2 == err.second.end()) {
                err.second.insert(qubit);
                err.first.erase(it1);
            }
        } else if (it2 != err.second.end()) {
            err.first.insert(qubit);
            err.second.erase(it2);
        }
    }
}
// Phase gate
// Adds an additional Z error for every X error
void FrameSimulator::s(int qubit)
{
    for (auto &[shot, err] : errors) {
        if (err.first.find(qubit) != err.first.end())
            flip_error(shot, qubit, ERROR_Z);
    }
}
// Sqrt(X) gate
// Adds an additional X error for every Z error
void FrameSimulator::sx(int qubit)
{
    for (auto &[shot, err] : errors) {
        if (err.second.find(qubit) != err.second.end())
            flip_error(shot, qubit, ERROR_X);
    }
}
// CNOT gate
// Propagates X errors from control to target, and Z errors from target to control
void FrameSimulator::cx(int control, int target)
{
    for (auto it = errors.begin(); it != errors.end(); ) {
        size_t shot = it->first;
        auto &err = it->second;
        if (err.first.find(control) != err.first.end())
            flip_error(it, target, ERROR_X);
        else
            ++it;
    }
    for (auto it = errors.begin(); it != errors.end(); ) {
        size_t shot = it->first;
        auto &err = it->second;
        if (err.second.find(target) != err.second.end())
            flip_error(it, control, ERROR_Z);
        else
            ++it;
    }
}
// CZ gate
void FrameSimulator::cz(int q1, int q2)
{
    for (auto it = errors.begin(); it != errors.end(); ) {
        size_t shot = it->first;
        auto &err = it->second;
        if (err.first.find(q1) != err.first.end())
            flip_error(it, q2, ERROR_Z);
        else
            ++it;
    }
    for (auto it = errors.begin(); it != errors.end(); ) {
        size_t shot = it->first;
        auto &err = it->second;
        if (err.first.find(q2) != err.first.end())
            flip_error(it, q1, ERROR_Z);
        else
            ++it;
    }
}
void FrameSimulator::sxx(int q1, int q2)
{
    for (auto &[shot, err] : errors) {
        if ((err.second.find(q1) != err.second.end()) ^ (err.second.find(q2) != err.second.end())) {
            flip_error(shot, q1, ERROR_X);
            flip_error(shot, q2, ERROR_X);
        }
    }
}
void FrameSimulator::szz(int q1, int q2)
{
    for (auto &[shot, err] : errors) {
        if ((err.first.find(q1) != err.first.end()) ^ (err.first.find(q2) != err.first.end())) {
            flip_error(shot, q1, ERROR_Z);
            flip_error(shot, q2, ERROR_Z);
        }
    }
}
// Measurement in the X basis
// Flipped if there is a Z error in the qubit
void FrameSimulator::mx(int qubit, MeasurementTag tag)
{
#ifdef RANDOMIZE_FLIPS
    for (size_t i=0; i<num_shots; i++) {
        if (randomizer(rng)) flip_error(i, qubit, ERROR_X);
    }
#endif
    for (auto &[shot, err] : errors) {
        if (err.second.find(qubit) != err.second.end())
            qubit_measurement_results[shot].results[qubit].insert(tag);
    }
}
// Measurement in the Z basis
// Flipped if there is a X error in the qubit
void FrameSimulator::mz(int qubit, MeasurementTag tag)
{
#ifdef RANDOMIZE_FLIPS
    for (size_t i=0; i<num_shots; i++) {
        if (randomizer(rng)) flip_error(i, qubit, ERROR_Z);
    }
#endif
    for (auto &[shot, err] : errors) {
        if (err.first.find(qubit) != err.first.end())
            qubit_measurement_results[shot].results[qubit].insert(tag);
    }
}
// Qubit preparation in |+>
// Resets Z errors in the qubit
void FrameSimulator::rx(int qubit)
{
#ifdef RANDOMIZE_FLIPS
    for (size_t i=0; i<num_shots; i++) {
        if (randomizer(rng)) flip_error(i, qubit, ERROR_X);
    }
#endif
    for (auto it = errors.begin(); it != errors.end(); ) {
        reset_error(it, qubit, ERROR_Z);
    }
}
// Qubit preparation in |0>
// Resets X errors in the qubit
void FrameSimulator::rz(int qubit)
{
#ifdef RANDOMIZE_FLIPS
    for (size_t i=0; i<num_shots; i++) {
        if (randomizer(rng)) flip_error(i, qubit, ERROR_Z);
    }
#endif
    for (auto it = errors.begin(); it != errors.end(); ) {
        reset_error(it, qubit, ERROR_X);
    }
}
// Runs a circuit node, which consists in a deterministic circuit
// and a function which determines which circuit goes afterwards
void FrameSimulator::run(std::shared_ptr<CircuitNode> node)
{
    //std::cout<<node->name<<std::endl;
    // Run the circuit
    BaseFrameSimulator::run(node->circuit);
    
    /*if (node->error_corrections && node->next_node_index)
        abort();*/
    // Apply error corrections for this round, if any
    if (node->error_corrections) {
        for (auto it = qubit_measurement_results.begin(); it != qubit_measurement_results.end(); ) {
            size_t shot = it->first;
            auto &res = it->second;
            auto corr = node->error_corrections(res);
            for (int errx : corr.first) {
                flip_error(shot, errx, ERROR_X);
            }
            for (int errz : corr.second) {
                flip_error(shot, errz, ERROR_Z);
            }
            if (res.results.empty())
                it = qubit_measurement_results.erase(it);
            else
                ++it;
        }
    }
    if (node->childs.size() <= 1 && !node->next_node_index) {
        if (!node->childs.empty())
            run(node->childs[0]);
        return;
    }
    // Sort shots by next circuit index
    std::map<int, std::vector<size_t>> branch_shots;
    // Process shots for which an error has been detected
    for (auto &[shot, res] : qubit_measurement_results) {
        int branch = 0;
        if (node->next_node_index)
            branch = node->next_node_index(res);
        branch_shots[branch].push_back(shot);
    }
    // Process noisy runs with undetected errors
    for (auto &[shot, err] : errors) {
        if (qubit_measurement_results.find(shot) == qubit_measurement_results.end()) {
            branch_shots[0].push_back(shot);
        }
    }
    int processed_shots=0;
    for (auto &kvp : branch_shots) {
        processed_shots += kvp.second.size();
    }
    //std::cout<<node->name<<" "<<errors.size()<<std::endl;
    //std::cout<<num_shots<<" "<<num_shots-processed_shots<<std::endl;
    auto old_err = std::move(errors);
    auto old_res = std::move(qubit_measurement_results);
    errors.clear();
    qubit_measurement_results.clear();
    int start = 0;
    // Run next circuits
    for (int i=0; i<node->childs.size() || i==0; i++) {
        auto &shots = branch_shots[i];
        int nshots = shots.size();
        // For the first branch, include shots which have no errors nor detected measurements
        if (i == 0)
            nshots += num_shots-processed_shots;
        if (nshots == 0)
            continue;
        if (i < node->childs.size() && node->childs[i]) {
            //if (i == 0 && num_shots-processed_shots > shots.size()) {
                FrameSimulator sim(nshots, rng);
                // Setup new simulation with initial states equal to end state of the recently run simulation
                // Fill with shots which require this branch
                for (int j=0; j<shots.size(); j++) {
                    auto it1 = old_err.find(shots[j]);
                    if (it1 != old_err.end())
                        sim.errors[j] = std::move(it1->second);
                    auto it2 = old_res.find(shots[j]);
                    if (it2 != old_res.end())
                        sim.qubit_measurement_results[j] = std::move(it2->second);    
                }
                sim.error = error;
                sim.current_tick = current_tick+1;
                // Run next circuit
                sim.run(node->childs[i]);
                nshots = sim.num_shots;
                // Move results and errors from the new simulation back to the original simulation
                for (auto &[index,err] : sim.errors) {
                    errors[index+start] = std::move(err);
                }
                for (auto &[index,res] : sim.qubit_measurement_results) {
                    qubit_measurement_results[index+start] = std::move(res);
                }
            /*} else {
                DenseFrameSimulator sim(0, rng);
                for (int j=0; j<shots.size(); j++) {
                    auto it1 = old_err.find(shots[j]);
                    if (it1 != old_err.end())
                        sim.errors[j] = std::move(it1->second);
                    auto it2 = old_res.find(shots[j]);
                    if (it2 != old_res.end())
                        sim.qubit_measurement_results[j] = std::move(it2->second);
                }
                for (int j=shots.size(); j<nshots; j++) {
                    sim.errors
                }
            }*/
        } else {
            // If there is no following circuit, just reorganize errors and sort them by branch
            for (int j=0; j<shots.size(); j++) {
                auto it1 = old_err.find(shots[j]);
                if (it1 != old_err.end())
                    errors[j+start] = std::move(it1->second);
                auto it2 = old_res.find(shots[j]);
                if (it2 != old_res.end())
                    qubit_measurement_results[j+start] = std::move(it2->second);
            }
        }
        start += nshots;
    }
    num_shots = start;
}