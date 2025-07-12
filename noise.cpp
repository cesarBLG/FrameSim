#include "noise.h"
#include <set>
std::set<InstructionType> inst1q({InstructionType::H, InstructionType::X, InstructionType::Y, InstructionType::Z, InstructionType::I, InstructionType::SX, InstructionType::SXDG, InstructionType::SY, InstructionType::SYDG, InstructionType::S, InstructionType::SDG});
std::set<InstructionType> inst2q({InstructionType::CX, InstructionType::CZ, InstructionType::SXX, InstructionType::SXXDG, InstructionType::SZZ, InstructionType::SZZDG});
std::set<InstructionType> instm({InstructionType::MX, InstructionType::MY, InstructionType::MZ, InstructionType::RX, InstructionType::RY, InstructionType::RZ});
std::set<InstructionType> instreset({InstructionType::RX, InstructionType::RY, InstructionType::RZ});
std::set<InstructionType> instnoise({InstructionType::DEPOLARIZE1, InstructionType::DEPOLARIZE2, InstructionType::DEPOLARIZE, InstructionType::X_ERROR, InstructionType::Y_ERROR, InstructionType::Z_ERROR});
Circuit DepolarizingModel::noisy_circuit(Circuit &circ)
{
    Circuit newcirc;
    std::set<int> gated_qubits;
    std::set<int> entangled_qubits;
    std::set<int> measured_qubits;
    for (auto it = circ.instructions.begin(); it != circ.instructions.end(); ++it) {
        Circuit pre;
        Circuit mid;
        Circuit post;
        auto &inst = *it;
        mid.append(inst);
        if (inst.type == InstructionType::TICK) {
            std::vector<int> idle1_qubits;
            std::vector<int> idle2_qubits;
            std::vector<int> idlem_qubits;
            if (!gated_qubits.empty()) {
                int start1 = 0;
                for (auto it2 = gated_qubits.begin(); it2 != gated_qubits.end(); ++it2) {
                    for (int i=start1; i<*it2; i++)
                        idle1_qubits.push_back(i);
                    start1 = *it2+1;
                }
                for (int i=start1; i<circ.num_qubits; i++)
                    idle1_qubits.push_back(i);
            }
            if (!entangled_qubits.empty()) {
                int start2 = 0;
                for (auto it2 = entangled_qubits.begin(); it2 != entangled_qubits.end(); ++it2) {
                    for (int i=start2; i<*it2; i++)
                        idle2_qubits.push_back(i);
                    start2 = *it2+1;
                }
                for (int i=start2; i<circ.num_qubits; i++)
                    idle2_qubits.push_back(i);
            }
            if (!measured_qubits.empty()/* && measured_qubits.size() > 5*/) { // FIXME
                int startm = 0;
                for (auto it2 = measured_qubits.begin(); it2 != measured_qubits.end(); ++it2) {
                    for (int i=startm; i<*it2; i++)
                        idlem_qubits.push_back(i);
                    startm = *it2+1;
                }
                for (int i=startm; i<circ.num_qubits; i++)
                    idlem_qubits.push_back(i);
            }
            gated_qubits.clear();
            entangled_qubits.clear();
            measured_qubits.clear();
            /*if (pgate > 0)
                pre.append(Instruction(InstructionType::DEPOLARIZE1, idle1_qubits, pidle));
            if (pidle-pgate > 0)
                pre.append(Instruction(InstructionType::DEPOLARIZE1, idle2_qubits, pidle-pgate));*/
            /*if (pgate > 0 && pidle > 0 && idle2_qubits.empty() && idlem_qubits.empty() && !idle1_qubits.empty())
                pre.append(Instruction(InstructionType::DEPOLARIZE1, idle1_qubits, std::min(pidle, pgate)));*/
            if (pidle > 0 && !idle2_qubits.empty())
                pre.append(Instruction(bias ? InstructionType::Z_ERROR : InstructionType::DEPOLARIZE1, idle2_qubits, pidle)); // Only decoherence
            if (deltapidle > 0 && !idlem_qubits.empty())
                pre.append(Instruction(InstructionType::DEPOLARIZE1, idlem_qubits, deltapidle));
        } else if (inst.type == InstructionType::MX || inst.type == InstructionType::MY || inst.type == InstructionType::MZ) {
            if (pm > 0)
                pre.append(Instruction(inst.type == InstructionType::MX ? InstructionType::Z_ERROR : InstructionType::X_ERROR, inst.targets, pm));
        } else if (inst.type == InstructionType::RX || inst.type == InstructionType::RY || inst.type == InstructionType::RZ) {
            if (pm > 0)
                post.append(Instruction(inst.type == InstructionType::RX ? InstructionType::Z_ERROR : InstructionType::X_ERROR, inst.targets, pm));
        } else if (inst.type == InstructionType::CX || inst.type == InstructionType::CY || inst.type == InstructionType::CZ) {
            if (pcnot > 0)
                post.append(Instruction(InstructionType::DEPOLARIZE2, inst.targets, pcnot));
        } else if (inst.type == InstructionType::SXX || inst.type == InstructionType::SXXDG || inst.type == InstructionType::SZZ || inst.type == InstructionType::SZZDG) {
            if (pcnot > 0)
                post.append(Instruction(InstructionType::DEPOLARIZE, inst.targets, pcnot));
        } else if (inst1q.find(inst.type) != inst1q.end()) {
            if (pgate > 0)
                post.append(Instruction(InstructionType::DEPOLARIZE1, inst.targets, pgate));
        }
        for (int q : inst.targets) {
            if (instnoise.find(inst.type) != instnoise.end())
                break;
            if (entangled_qubits.find(q) != entangled_qubits.end() || measured_qubits.find(q) != measured_qubits.end()) {
                abort();
            }
            if (instm.find(inst.type) != instm.end()) {
                measured_qubits.insert(q);
                entangled_qubits.insert(q);
                gated_qubits.insert(q);
            } else if (inst2q.find(inst.type) != inst2q.end()) {
                entangled_qubits.insert(q);
                gated_qubits.insert(q);
            } else if (inst1q.find(inst.type) != inst1q.end()) {
                gated_qubits.insert(q);
            }
        }
        newcirc += pre + mid + post;
    }
    return std::move(newcirc);
}
Circuit GeneralDepolarizingModel::noisy_circuit(Circuit &circ)
{
    Circuit newcirc;
    std::set<int> gated_qubits;
    std::set<int> entangled_qubits;
    std::set<int> measured_qubits;
    std::map<int,double> used_time;
    double cooling_time = 0;
    for (auto it = circ.instructions.begin(); it != circ.instructions.end(); ++it) {
        Circuit pre;
        Circuit mid;
        Circuit post;
        auto &inst = *it;
        mid.append(inst);
        if (inst.type == InstructionType::TICK) {
            double maxtime = 0;
            for (auto &kvp : used_time) {
                if (kvp.second > maxtime)
                    maxtime = kvp.second;
            }
            maxtime += cooling_time;
            for (int i=0; i<circ.num_qubits; i++) {
                double time = maxtime-used_time[i];
                if (time <= 0 || (T1 == 0 && T2 == 0)) {
                    continue;
                } else if (T1 == 0) {
                    pre.append(Instruction(InstructionType::Z_ERROR, {i}, time/2/T2));
                } else if (T1 == T2) {
                    pre.append(Instruction(InstructionType::DEPOLARIZE1, {i}, time/2*(1/(2*T1)+1/T2)));
                } else {
                    double px = time/4/T1;
                    double pz = time/2*(1/T2-0.5/T1);
                    pre.append(Instruction(InstructionType::PAULI1, {i}, {px, px, pz}));
                }
            }
            used_time.clear();
            gated_qubits.clear();
            entangled_qubits.clear();
            measured_qubits.clear();
            cooling_time = 0;
        } else {
            double p = errors[inst.type];
            if (p > 0) {
                if (inst1q.find(inst.type) != inst1q.end()) {
                    post.append(Instruction(InstructionType::DEPOLARIZE1, inst.targets, p));
                } else if (inst.type == InstructionType::MX || inst.type == InstructionType::MY || inst.type == InstructionType::MZ) {
                    pre.append(Instruction(inst.type == InstructionType::MX ? InstructionType::Z_ERROR : InstructionType::X_ERROR, inst.targets, p));
                } else if (inst.type == InstructionType::RX || inst.type == InstructionType::RY || inst.type == InstructionType::RZ) {
                    post.append(Instruction(inst.type == InstructionType::RX ? InstructionType::Z_ERROR : InstructionType::X_ERROR, inst.targets, p));
                } else if (inst.type == InstructionType::CX || inst.type == InstructionType::CY || inst.type == InstructionType::CZ) {
                    post.append(Instruction(InstructionType::DEPOLARIZE2, inst.targets, p));
                } else if (inst.type == InstructionType::SXX || inst.type == InstructionType::SXXDG || inst.type == InstructionType::SZZ || inst.type == InstructionType::SZZDG) {
                    post.append(Instruction(InstructionType::DEPOLARIZE, inst.targets, p));
                }
            }
            if (inst.type == InstructionType::DELAY) {
                auto it = used_time.find(-1);
                double delay = delay_times[*inst.label];
                if (it == used_time.end() || it->second < delay)
                    used_time[-1] = delay;
                if (delay_cooling_times.find(*inst.label) != delay_cooling_times.end()) {
                    double cool = delay_cooling_times[*inst.label];
                    if (cool > cooling_time)
                        cooling_time = cool;
                }
            } else {
                if (cooling_times.find(inst.type) != cooling_times.end()) {
                    double cool = cooling_times[inst.type];
                    if (cool > cooling_time)
                        cooling_time = cool;
                }
            }
            for (int q : inst.targets) {
                if (instnoise.find(inst.type) != instnoise.end())
                    break;
                used_time[q] += times[inst.type];
                if (entangled_qubits.find(q) != entangled_qubits.end() || measured_qubits.find(q) != measured_qubits.end())
                    abort();
                else if (instm.find(inst.type) != instm.end()) {
                    measured_qubits.insert(q);
                    entangled_qubits.insert(q);
                    gated_qubits.insert(q);
                } else if (inst2q.find(inst.type) != inst2q.end()) {
                    entangled_qubits.insert(q);
                    gated_qubits.insert(q);
                } else if (inst1q.find(inst.type) != inst1q.end()) {
                    gated_qubits.insert(q);
                }
            }
        }
        newcirc += pre + mid + post;
    }
    return std::move(newcirc);
}
Circuit MidCircuitNoiseModel::noisy_circuit(Circuit &circ)
{
    Circuit newcirc;
    std::set<int> gated_qubits;
    std::set<int> entangled_qubits;
    std::set<int> measured_qubits;
    for (auto it = circ.instructions.begin(); it != circ.instructions.end(); ++it) {
        Circuit pre;
        Circuit mid;
        Circuit post;
        auto &inst = *it;
        mid.append(inst);
        if (inst.type == InstructionType::TICK) {
            std::vector<int> idle1_qubits;
            std::vector<int> idle2_qubits;
            std::vector<int> idlem_qubits;
            if (!gated_qubits.empty()) {
                int start1 = 0;
                for (auto it2 = gated_qubits.begin(); it2 != gated_qubits.end(); ++it2) {
                    for (int i=start1; i<*it2; i++)
                        idle1_qubits.push_back(i);
                    start1 = *it2+1;
                }
                for (int i=start1; i<circ.num_qubits; i++)
                    idle1_qubits.push_back(i);
            }
            if (!entangled_qubits.empty()) {
                int start2 = 0;
                for (auto it2 = entangled_qubits.begin(); it2 != entangled_qubits.end(); ++it2) {
                    for (int i=start2; i<*it2; i++)
                        idle2_qubits.push_back(i);
                    start2 = *it2+1;
                }
                for (int i=start2; i<circ.num_qubits; i++)
                    idle2_qubits.push_back(i);
            }
            if (!measured_qubits.empty()) {
                int startm = 0;
                for (auto it2 = measured_qubits.begin(); it2 != measured_qubits.end(); ++it2) {
                    for (int i=startm; i<*it2; i++)
                        idlem_qubits.push_back(i);
                    startm = *it2+1;
                }
                for (int i=startm; i<circ.num_qubits; i++)
                    idlem_qubits.push_back(i);
            }
            if (!idlem_qubits.empty()) {
                if (measured_qubits.size() > 1)
                    pre.append(err_midcirc.get_instruction(idlem_qubits));
                //else
                    //pre.append(Error::DelayError(t_1q, 0, T2).get_instruction(idlem_qubits));
            } else {
                if (!idle1_qubits.empty())
                    pre.append(Error::DelayError(t_1q, 0, T2).get_instruction(idle1_qubits));
                if (!idle2_qubits.empty())
                    pre.append(Error::DelayError(t_2q-t_1q, 0, T2).get_instruction(idle2_qubits));
            }
            gated_qubits.clear();
            entangled_qubits.clear();
            measured_qubits.clear();
        } else if (inst.type == InstructionType::MX || inst.type == InstructionType::MY || inst.type == InstructionType::MZ) {
            if (err_m > 0)
                pre.append(Instruction(inst.type == InstructionType::MX ? InstructionType::Z_ERROR : InstructionType::X_ERROR, inst.targets, err_m));
        } else if (inst.type == InstructionType::RX || inst.type == InstructionType::RY || inst.type == InstructionType::RZ) {
            if (err_m > 0)
                post.append(Instruction(inst.type == InstructionType::RX ? InstructionType::Z_ERROR : InstructionType::X_ERROR, inst.targets, err_m));
        } else if (inst.type == InstructionType::CX || inst.type == InstructionType::CY || inst.type == InstructionType::CZ) {
            if (err_2q > 0)
                post.append(Instruction(InstructionType::DEPOLARIZE2, inst.targets, err_2q));
        } else if (inst.type == InstructionType::SXX || inst.type == InstructionType::SXXDG || inst.type == InstructionType::SZZ || inst.type == InstructionType::SZZDG) {
            if (err_2q > 0)
                post.append(Instruction(InstructionType::DEPOLARIZE, inst.targets, err_2q));
        } else if (inst1q.find(inst.type) != inst1q.end()) {
            if (err_1q > 0)
                post.append(Instruction(InstructionType::DEPOLARIZE1, inst.targets, err_1q));
        }
        for (int q : inst.targets) {
            if (instnoise.find(inst.type) != instnoise.end())
                break;
            if (entangled_qubits.find(q) != entangled_qubits.end() || measured_qubits.find(q) != measured_qubits.end())
                abort();
            if (instm.find(inst.type) != instm.end()) {
                if (instreset.find(inst.type) == instreset.end())
                    measured_qubits.insert(q);
                //entangled_qubits.insert(q);
                gated_qubits.insert(q);
            } else if (inst2q.find(inst.type) != inst2q.end()) {
                entangled_qubits.insert(q);
                gated_qubits.insert(q);
            } else if (inst1q.find(inst.type) != inst1q.end()) {
                gated_qubits.insert(q);
            }
        }
        newcirc += pre + mid + post;
    }
    return std::move(newcirc);
}
void apply_noise_to_nodes(std::shared_ptr<CircuitNode> node0, NoiseModel &noise, std::set<std::shared_ptr<CircuitNode>> &visited)
{
    if (visited.find(node0) != visited.end())
        return;
    visited.insert(node0);
    node0->circuit = noise.noisy_circuit(node0->circuit);
    for (auto node : node0->childs) {
        if (node != nullptr)
            apply_noise_to_nodes(node, noise, visited);
    }
}
// Applies the specified noise model to all circuits in a tree
void apply_noise_to_nodes(std::shared_ptr<CircuitNode> node0, NoiseModel &noise)
{
    std::set<std::shared_ptr<CircuitNode>> visited;
    apply_noise_to_nodes(node0, noise, visited);
}
// Example depolarizing noise model for ion traps, with current and anticipated error rates
// TODO: move it to an examples folder, this is not part of the simulator itself
GeneralDepolarizingModel::GeneralDepolarizingModel(double alpha, double T2) : T1(0), T2(T2)
{
    double beta = 1-alpha;
    for (auto &instr : inst1q)
    {
        errors[instr] = 0.0036*alpha+1e-5*beta;
        times[instr] = 25e-6*alpha+1e-6*beta;
    }
    for (auto &instr : inst2q)
    {
        errors[instr] = 0.027*alpha + 2e-4*beta;
        times[instr] = 322e-6*alpha + 15e-6*beta;
    }
    for (auto &instr : instm)
    {
        errors[instr] = 3e-3*alpha + 1e-4*beta;
        times[instr] = 400e-6*alpha+30e-6*beta;
        if (instreset.find(instr) == instreset.end())
            cooling_times[instr] = 100e-6*alpha+50e-6*alpha;
    }
    for (auto &instr : instreset)
    {
        errors[instr] = 0.003;
        times[instr] = 50e-6*alpha+10e-6*beta;
    }
    delay_times["SPLIT"] = delay_times["MERGE"] = 80e-6*alpha+30e-6*beta;
    delay_times["JUNCTION_TRANSPORT"] = 200e-6*alpha+100e-6*beta;
    delay_times["ROTATION"] = 150e-6*alpha+20e-6*beta;
    delay_cooling_times["JUNCTION_TRANSPORT"] = 100e-6*alpha+25e-6*beta;
    delay_cooling_times["ROTATION"] = 100e-6*alpha+25e-6*beta;
    //delays["RECOOLING"] = 400e-6*alpha+100e-6*beta;
}