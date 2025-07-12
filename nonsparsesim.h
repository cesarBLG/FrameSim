#pragma once
#include <vector>
#include <cstdint>
#include <map>
#include <set>
#include "simulator_base.h"
// Stores whether a flip is recorded for each shot
// Stored as a uint64_t array where each bit corresponds to a different shot
struct ErrorTable
{
    size_t nshots=0;
    std::vector<uint64_t> shots;
    ErrorTable& operator^=(const ErrorTable &o)
    {
        for (int i=0; i<shots.size(); i++) {
            shots[i] ^= o.shots[i];
        }
        return *this;
    }
    ErrorTable operator^(const ErrorTable &o) const
    {
        ErrorTable t = *this;
        for (int i=0; i<shots.size(); i++) {
            t.shots[i] ^= o.shots[i];
        }
        return t;
    }
    void flip(size_t shot)
    {
        shots[shot>>6] ^= UINT64_C(1)<<(shot&63);
    }
    void reset()
    {
        for (int i=0; i<shots.size(); i++) {
            shots[i] = 0;
        }
    }
    bool flipped(size_t shot)
    {
        return (shots[shot>>6] & (UINT64_C(1)<<(shot&63))) != 0;
    }
    bool reset_flipped(size_t shot, bool val=false)
    {
        uint64_t maj = shot>>6;
        uint64_t min = shot & 63;
        bool was_flipped = (shots[maj]&(UINT64_C(1)<<min)) != 0;
        shots[maj] ^= ((uint64_t)(val^was_flipped))<<min;
        return was_flipped;
    }
    void append(ErrorTable &&t)
    {
        while ((nshots & 63) && t.nshots) {
            reset_flipped(nshots++, t.flipped(--t.nshots));
            if ((t.nshots & 63) == 0)
                t.shots.pop_back();
        }
        shots.insert(shots.end(), t.shots.begin(), t.shots.end());
        nshots += t.nshots;
    }
    void append(bool err)
    {
        uint64_t maj = nshots>>6;
        uint64_t min = nshots&63;
        if (min == 0)
            shots.push_back(0);
        shots[maj] ^= (shots[maj]&(UINT64_C(1)<<min))|(((uint64_t)err)<<min);
        ++nshots;
    }
};
// Stores an error table for every measurement, indicating whether it has been flipped
class MeasurementResultsDense : public MeasurementResults
{
    std::map<int, std::map<MeasurementTag, ErrorTable>> &qubit_measurement_results;
    size_t shot;
    size_t num_shots;
    public:
    MeasurementResultsDense(std::map<int, std::map<MeasurementTag, ErrorTable>> &results, size_t shot, size_t nshots) : qubit_measurement_results(results), shot(shot), num_shots(nshots) {}
    bool is_flipped(int qubit, const MeasurementTag &tag) override
    {
        auto &tab = qubit_measurement_results[qubit][tag];
        tab.shots.resize(((num_shots-1)>>6)+1);
        return tab.flipped(shot);
    }
    bool reset_flipped(int qubit, const MeasurementTag &tag) override
    {
        auto &tab = qubit_measurement_results[qubit][tag];
        tab.shots.resize(((num_shots-1)>>6)+1);
        return tab.reset_flipped(shot);
    }
    void flip(int qubit, const MeasurementTag &tag) override
    {
        auto &tab = qubit_measurement_results[qubit][tag];
        tab.shots.resize(((num_shots-1)>>6)+1);
        tab.flip(shot);
    }
};
// Frame simulator using a dense flip representation
// Works similarly to a non-optimized version of Stim, but with support
// for dynamic circuits. It is recommended to use FrameSimulator instead,
// unless shot randomization is required to check the correctness of a protocol
class DenseFrameSimulator : public BaseFrameSimulator
{
    public:
    std::map<int,std::pair<ErrorTable,ErrorTable>> errors;
    // Tracks measurement results which have been flipped due to an error
    std::map<int, std::map<MeasurementTag, ErrorTable>> qubit_measurement_results;
    DenseFrameSimulator(size_t num_shots, std::mt19937_64 &rng) : BaseFrameSimulator(num_shots, rng) {}
    void h(int qubit) override;
    void s(int qubit) override;
    void sx(int qubit) override;
    void cx(int control, int target) override;
    void cz(int q1, int q2) override;
    void sxx(int q1, int q2) override;
    void szz(int q1, int q2) override;
    void mx(int qubit, MeasurementTag tag) override;
    void mz(int qubit, MeasurementTag tag) override;
    void rx(int qubit) override;
    void rz(int qubit) override;
    void run(Circuit &circ) override;
    void run(std::shared_ptr<CircuitNode> node) override;
    void flip_error(size_t shot, int qubit, int type) override;
};