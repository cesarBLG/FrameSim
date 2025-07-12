#pragma once
#include <random>
#include "circuit.h"
#define ERROR_X 1
#define ERROR_Z 2
// Base class for a frame simulator which is capable to run circuit trees given the initial node
class BaseFrameSimulator
{
    protected:
    std::mt19937_64 &rng;
    // Randomly determines the index of next faulty shot
    std::bernoulli_distribution randomizer;
    // Randomly picks a number between 1 and 3 for single qubit depolarizing
    std::uniform_int_distribution<> depolarizer1;
    // Randomly picks a number between 1 and 15 for two-qubit depolarizing
    std::uniform_int_distribution<> depolarizer2;
    // Current timestep of the simulation
    int current_tick=0;
    // Number of shots
    size_t num_shots;
    bool error=false;
    public:
    BaseFrameSimulator(size_t num_shots, std::mt19937_64 &rng) : num_shots(num_shots), rng(rng), depolarizer1(1,3), depolarizer2(1,15) {}
    virtual ~BaseFrameSimulator() = default;
    virtual void h(int qubit)=0;
    virtual void s(int qubit)=0;
    virtual void sx(int qubit)=0;
    virtual void sxx(int q1, int q2)=0;
    virtual void szz(int q1, int q2)=0;
    virtual void cx(int control, int target)=0;
    virtual void cz(int q1, int q2)=0;
    virtual void mx(int qubit, MeasurementTag tag)=0;
    virtual void mz(int qubit, MeasurementTag tag)=0;
    virtual void rx(int qubit)=0;
    virtual void rz(int qubit)=0;
    virtual void x_error(int qubit, double p);
    virtual void y_error(int qubit, double p);
    virtual void z_error(int qubit, double p);
    virtual void depolarize(std::vector<int> &qubits, double p);
    virtual void depolarize1(int qubit, double p);
    virtual void depolarize2(int control, int target, double p);
    virtual void pauli1(int qubit, std::vector<double> &p);
    virtual void pauli2(int q1, int q2, std::vector<double> &p);
    virtual void run(Instruction &instruction);
    virtual void run(Circuit &circ);
    virtual void run(std::shared_ptr<CircuitNode> node)=0;
    virtual void flip_error(size_t shot, int qubit, int type)=0;
    inline size_t get_num_shots() { return num_shots; }
};