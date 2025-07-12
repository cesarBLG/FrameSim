#include "nonsparsesim.h"
#define RANDOMIZE_FLIPS
void DenseFrameSimulator::flip_error(size_t shot, int qubit, int type)
{
    auto &err = errors[qubit];
    if (type & ERROR_X) err.first.flip(shot);
    if (type & ERROR_Z) err.second.flip(shot);
}
// Hadamard gate
// Exchanges X and Z errors in the target qubit
void DenseFrameSimulator::h(int qubit)
{
    auto &errs = errors[qubit];
    auto tmp = std::move(errs.first);
    errs.first = std::move(errs.second);
    errs.second = std::move(tmp);
}
// Phase gate
// Adds an additional Z error for every X error
void DenseFrameSimulator::s(int qubit)
{
    auto &errs = errors[qubit];
    errs.second ^= errs.first;
}
// Phase gate
// Adds an additional X error for every Z error
void DenseFrameSimulator::sx(int qubit)
{
    auto &errs = errors[qubit];
    errs.first ^= errs.second;
}
// CNOT gate
// Propagates X errors from control to target, and Z errors from target to control
void DenseFrameSimulator::cx(int control, int target)
{
    auto &errs1 = errors[control];
    auto &errs2 = errors[target];
    errs2.first ^= errs1.first;
    errs1.second ^= errs2.second;
}
// CZ gate
void DenseFrameSimulator::cz(int q1, int q2)
{
    auto &errs1 = errors[q1];
    auto &errs2 = errors[q2];
    errs2.second ^= errs1.first;
    errs1.second ^= errs2.first;
}
// Rxx(pi/2) = Sqrt(XX) gate
void DenseFrameSimulator::sxx(int q1, int q2)
{
    auto &errs1 = errors[q1];
    auto &errs2 = errors[q2];
    auto tmp = errs1.second ^ errs2.second;
    errs1.first ^= tmp;
    errs2.first ^= tmp;
}
// Rzz(pi/2) = Sqrt(ZZ) gate
void DenseFrameSimulator::szz(int q1, int q2)
{
    auto &errs1 = errors[q1];
    auto &errs2 = errors[q2];
    auto tmp = errs1.first ^ errs2.first;
    errs1.second ^= tmp;
    errs2.second ^= tmp;
}
// Measurement in the X basis
// Flipped if there is a Z error in the qubit
void DenseFrameSimulator::mx(int qubit, MeasurementTag tag)
{
    auto &errs = errors[qubit];
#ifdef RANDOMIZE_FLIPS
    for (int i=0; i<errs.first.shots.size(); i++) {
        errs.first.shots[i] = rng();
    }
#endif
    qubit_measurement_results[qubit][tag] = errs.second;
}
// Measurement in the Z basis
// Flipped if there is a X error in the qubit
void DenseFrameSimulator::mz(int qubit, MeasurementTag tag)
{
    auto &errs = errors[qubit];
#ifdef RANDOMIZE_FLIPS
    for (int i=0; i<errs.second.shots.size(); i++) {
        errs.second.shots[i] = rng();
    }
#endif
    qubit_measurement_results[qubit][tag] = errs.first;
}
// Qubit preparation in |+>
// Resets Z errors in the qubit
void DenseFrameSimulator::rx(int qubit)
{
    auto &errs = errors[qubit];
#ifdef RANDOMIZE_FLIPS
    for (int i=0; i<errs.first.shots.size(); i++) {
        errs.first.shots[i] = rng();
    }
#endif
    errs.second.reset();
}
// Qubit preparation in |0>
// Resets X errors in the qubit
void DenseFrameSimulator::rz(int qubit)
{
    auto &errs = errors[qubit];
#ifdef RANDOMIZE_FLIPS
    for (int i=0; i<errs.second.shots.size(); i++) {
        errs.second.shots[i] = rng();
    }
#endif
    errs.first.reset();
}
// Runs a deterministic circuit
void DenseFrameSimulator::run(Circuit &circuit)
{
    uint64_t sz = ((num_shots-1)>>6)+1;
    for (int i=0; i<circuit.num_qubits; i++) {
        auto &err = errors[i];
        err.first.shots.resize(sz);
        err.first.nshots = num_shots;
        if (err.second.shots.size() < sz) {
            int old = err.second.shots.size();
            err.second.shots.resize(sz);
            err.second.nshots = num_shots;
#ifdef RANDOMIZE_FLIPS
            for (int j=old; j<err.second.shots.size(); j++) {
                err.second.shots[j] = rng();
            }
#endif
        }
    }
    BaseFrameSimulator::run(circuit);
}
// Recursively runs a circuit node, which consists in a deterministic circuit
// and a function which determines which circuit goes afterwards
void DenseFrameSimulator::run(std::shared_ptr<CircuitNode> node)
{
    // Run the circuit
    run(node->circuit);

    // Apply error corrections for this round, if any
    if (node->error_corrections) {
        for (size_t i=0; i<num_shots; i++) {
            auto res = MeasurementResultsDense(qubit_measurement_results, i, num_shots);
            auto corr = node->error_corrections(res);
            for (int errx : corr.first) {
                errors[errx].first.flip(i);
            }
            for (int errz : corr.second) {
                errors[errz].second.flip(i);
            }
        }
    }
    // If only one node, directly run it
    if (node->childs.size() <= 1 && !node->next_node_index) {
        if (!node->childs.empty())
            run(node->childs[0]);
        return;
    }
    // Go over all shots and classify them depending on the measurement outcomes
    auto old_err = std::move(errors);
    auto old_res = std::move(qubit_measurement_results);
    errors.clear();
    qubit_measurement_results.clear();
    std::map<int, DenseFrameSimulator*> sims;
    for (size_t i=0; i<num_shots; i++) {
        auto res = MeasurementResultsDense(old_res, i, num_shots);
        // Get the index of the next circuit for this shot
        int branch = 0;
        if (node->next_node_index)
            branch = node->next_node_index(res);
        if (branch < 0) // Shot discarded by post-selection
            continue;
        // Build a new simulator for every circuit
        auto it = sims.find(branch);
        if (it == sims.end()) {
            it = sims.insert({branch, new DenseFrameSimulator(0, rng)}).first;
        }
        // Copy to the relevant simulator instance the shot data
        auto *sim = it->second;
        for (auto &[qub, tab] : old_err) {
            sim->errors[qub].first.append(tab.first.flipped(i));
            sim->errors[qub].second.append(tab.second.flipped(i));
        }
        for (auto &[qub, res] : old_res) {
            for (auto &[tag, tab] : res) {
                sim->qubit_measurement_results[qub][tag].append(tab.flipped(i));
            }
        }
        sim->num_shots++;
    }
    num_shots = 0;
    // Run the next circuit for each shot
    for (auto [branch, sim] : sims) {
        sim->error = error;
        sim->current_tick = current_tick;
        if (node->childs[branch])
            sim->run(node->childs[branch]);
        
        // Copy back the shot data to this instance
        for (auto &[qub, tab] : sim->errors) {
            errors[qub].first.append(std::move(tab.first));
            errors[qub].second.append(std::move(tab.second));
        }
        for (auto &[qub, res] : sim->qubit_measurement_results) {
            for (auto &[tag, tab] : res) {
                qubit_measurement_results[qub][tag].append(std::move(tab));
            }
        }
        num_shots += sim->num_shots;
        delete sim;
    }
}