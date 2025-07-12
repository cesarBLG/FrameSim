#pragma once
#include "simulator.h"
#include <list>
// Base noisy model class. Contains a noisy_circuit() method that takes one circuit
// and  returns its noisy version according to the defined noise model
class NoiseModel
{
    public:
    double pcrosstalk=0;
    double pmidcirc=0;
    virtual Circuit noisy_circuit(Circuit &circ)
    {
        return circ;
    }
    virtual Circuit cx_with_crosstalk(int control, int target, std::set<int> neighbours)
    {
        Circuit circ;
        circ.append(InstructionType::CX, control, target);
        if (pcrosstalk == 0)
            return circ;
        for (int q : neighbours) {
            if (q == control || q == target)
                continue;
            circ.append(Instruction(InstructionType::DEPOLARIZE2, {control, q}, pcrosstalk));
            circ.append(Instruction(InstructionType::DEPOLARIZE2, {target, q}, pcrosstalk));
        }
        return circ;
    }
};
// Standard circuit level noise with different depolarizing rates for each operation
class DepolarizingModel : public NoiseModel
{
    double pidle;
    double deltapidle;
    double pgate;
    double pcnot;
    double pm;
    bool bias=false;
    public:
    DepolarizingModel(double p) : pidle(p), deltapidle(0), pgate(p), pcnot(p), pm(p)
    {

    }
    DepolarizingModel(double p, double q) : pidle(p), deltapidle(0), pgate(p), pcnot(q), pm(p)
    {

    }
    DepolarizingModel(double pidle, double deltapidle, double pgate, double pcnot, double pm) : pidle(pidle), deltapidle(deltapidle), pgate(pgate), pcnot(pcnot), pm(pm)
    {
        
    }
    Circuit noisy_circuit(Circuit &circ) override;
};
// General error channel definition, that provides an error instruction to
// the specified qubits according to the channel
struct Error
{
    InstructionType type;
    std::vector<double> error_rates;
    Error() = default;
    Error(InstructionType type, std::vector<double> error_rates) : type(type), error_rates(error_rates) {}
    Error(InstructionType type, double p) : type(type), error_rates({p}) {}
    Instruction get_instruction(std::vector<int> qubits)
    {
        return Instruction(type, qubits, error_rates);
    }
    static Error DelayError(double time, double T1, double T2)
    {
        if ((T1 > 0 && T2 > 2*T1) || T2 == 0) abort();
        if (T1 == 0) {
            return Error(InstructionType::Z_ERROR, time/2/T2);
        } else if (T1 == T2) { 
            return Error(InstructionType::DEPOLARIZE1, time/2*(1/(2*T1)+1/T2));
        } else {
            double px = time/4/T1;
            double pz = time/2*(1/T2-0.5/T1);
            return Error(InstructionType::PAULI1, {px, px, pz});
        }
    }
};
// Biased depolarizing noise model, with a special handling of idle errors during
// mid-circuit measurements
class MidCircuitNoiseModel : public NoiseModel
{
    double T2;
    Error err_midcirc;
    double err_m;
    double err_2q;
    double err_1q;
    double t_2q; 
    double t_1q;
    public:
    MidCircuitNoiseModel(double T2, double err_1q, double t_1q, double err_2q, double t_2q, double err_m, std::vector<double> err_midcirc) : 
        T2(T2), err_m(err_m), err_2q(err_2q), err_1q(err_1q), t_2q(t_2q), t_1q(t_1q)
    {
        this->err_midcirc = Error(InstructionType::PAULI1, err_midcirc);
    }
    Circuit noisy_circuit(Circuit &circ) override;
};
// Generic noise model with custom depolarizing instructions for each gate
// Timings of operations, as well as required cooling times can be provided
// to generate idle errors
class GeneralDepolarizingModel : public NoiseModel
{
    double T1;
    double T2;
    std::map<InstructionType,double> errors;
    std::map<InstructionType,double> times;
    //std::map<std::string,double> delay_errors;
    std::map<std::string,double> delay_times;
    std::map<InstructionType,double> cooling_times;
    std::map<std::string,double> delay_cooling_times;
    public:
    GeneralDepolarizingModel()
    {

    }
    GeneralDepolarizingModel(double T1, double T2, std::map<InstructionType,double> errors, std::map<InstructionType,double> times/*, std::map<std::string, double> delay_errors*/, std::map<std::string, double> delay_times, std::map<InstructionType,double> cooling_times, std::map<std::string, double> delay_cooling_times) : T1(T1), T2(T2), errors(errors), times(times), /*delay_errors(delay_errors), */delay_times(delay_times), cooling_times(cooling_times), delay_cooling_times(delay_cooling_times)
    {
        if (T1 > 0 && T2 > 2*T1) abort();
    }
    GeneralDepolarizingModel(double alpha, double T2=50e-3);
    Circuit noisy_circuit(Circuit &circ) override;
};
void apply_noise_to_nodes(std::shared_ptr<CircuitNode> node0, NoiseModel &noise, std::set<std::shared_ptr<CircuitNode>> &visited);
void apply_noise_to_nodes(std::shared_ptr<CircuitNode> node0, NoiseModel &noise);