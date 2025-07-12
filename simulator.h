#pragma once
#include <vector>
#include <cstdint>
#include <map>
#include <set>
#include "simulator_base.h"
struct MeasurementResultsSparse : MeasurementResults
{
    std::map<int,std::set<MeasurementTag>> results;
    bool is_flipped(int qubit, const MeasurementTag &tag) override
    {
        auto it = results.find(qubit);
        return it != results.end() && it->second.find(tag) != it->second.end();
    }
    bool reset_flipped(int qubit, const MeasurementTag &tag) override
    {
        auto it = results.find(qubit);
        if (it != results.end()) {
            bool exists = it->second.erase(tag);
            if (exists) {
                if (it->second.empty())
                    results.erase(it);
                return true;
            }
        }
        return false;
    }
    void flip(int qubit, const MeasurementTag &tag) override
    {
        auto it = results.find(qubit);
        if (it != results.end()) {
            bool exists = it->second.erase(tag);
            if (exists) {
                if (it->second.empty())
                    results.erase(it);
                return;
            }
        }
        results[qubit].insert(tag);
    }
};
// Pauli frame simulator: it propagates Pauli errors to the end of the circuit
// QEC circuits are in the Clifford group: Pauli errors are mapped to Pauli errors after each gate
// Much more efficient than a full stabilizer or statevector simulation
// Executes shots in paralell: all shots are updated after each instruction
// Uses a sparse representation: only shots with errors are stored
// Warning: in order to use a sparse shot approach, 
// error randomization on measured and reset qubits in the opposite basis is disabled
// It must be ensured that the (noiseless) circuit is correct by other means before it is simulated
// To enable randomization to check circuit validity, use #define RANDOMIZE_FLIPS in simulator.cpp
// or use the DenseFrameSimulator
class FrameSimulator : public BaseFrameSimulator
{
    // Tracks errors for each shot
    // Key: shot number
    // Value: list of qubits with X (first) and Z (second) errors
    std::map<size_t, std::pair<std::set<int>,std::set<int>>> errors;
    public:
    // Tracks measurement results which have been flipped due to an error
    // Key: shot number
    // Value: list of flipped measurements, ordered by qubit
    std::map<size_t, MeasurementResultsSparse> qubit_measurement_results;
    FrameSimulator(size_t num_shots, std::mt19937_64 &rng) : BaseFrameSimulator(num_shots, rng) {}
    void h(int qubit) override;
    void s(int qubit) override;
    void sx(int qubit) override;
    void sxx(int q1, int q2) override;
    void szz(int q1, int q2) override;
    void cx(int control, int target) override;
    void cz(int q1, int q2) override;
    void mx(int qubit, MeasurementTag tag) override;
    void mz(int qubit, MeasurementTag tag) override;
    void rx(int qubit) override;
    void rz(int qubit) override;
    void run(std::shared_ptr<CircuitNode> node) override;
    void reset_error(std::map<size_t,std::pair<std::set<int>,std::set<int>>>::iterator &it, int qubit, int type);
    void flip_error(std::map<size_t,std::pair<std::set<int>,std::set<int>>>::iterator &it, int qubit, int type);
    void flip_error(size_t shot, int qubit, int type) override;
    inline size_t get_num_shots() { return num_shots; }
};