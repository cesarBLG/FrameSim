#include "circuit.h"
#include <string>
#include <iostream>
// String representation of the InstructionType enum
std::string InstructionNames[] = {
    "I",
    "X",
    "Y",
    "Z",
    "H",
    "S",
    "SDG",
    "SX",
    "SXDG",
    "SY",
    "SYDG",
    "CX",
    "CY",
    "CZ",
    "SXX",
    "SXXDG",
    "SZZ",
    "SZZDG",
    "MX",
    "MY",
    "MZ",
    "RX",
    "RY",
    "RZ",
    "DEPOLARIZE",
    "DEPOLARIZE1",
    "DEPOLARIZE2",
    "X_ERROR",
    "Y_ERROR",
    "Z_ERROR",
    "PAULI1",
    "PAULI2",
    "DELAY",
    "TICK",
};
// Prints a circuit instruction to the stream
std::ostream &operator<<(std::ostream &os, const Instruction &instr)
{
    os << InstructionNames[(int)instr.type];
    if (!instr.p.empty()) {
        os << '(';
        os << instr.p[0];
        for (int i=1; i<instr.p.size(); i++) {
            os << ',';
            os << instr.p[i];
        }
        os << ')';
    }
    if (instr.label)
        os << ' ' << *instr.label;
    for (int t : instr.targets) {
        os << ' ' << t;
    }
    return os;
}
// Prints a circuit to the stream
std::ostream &operator<<(std::ostream &os, const Circuit &circ)
{
    for (auto &instr : circ.instructions) {
        os << instr << '\n';
    }
    return os;
}
// Combines the operations from two circuits, keeping instructions in
// the original timestep, as specified by TICK instructions
Circuit merge_circuits(Circuit c1, Circuit c2)
{
    Circuit c;
    auto it1 = c1.instructions.begin();
    auto it2 = c2.instructions.begin();
    while (it1 != c1.instructions.end() || it2 != c2.instructions.end()) {
        while (it1 != c1.instructions.end()) {
            if (it1->type == InstructionType::TICK) {
                ++it1;
                break;
            } else {
                c.append(*it1);
                ++it1;
            }
        }
        while (it2 != c2.instructions.end()) {
            if (it2->type == InstructionType::TICK) {
                ++it2;
                break;
            } else {
                c.append(*it2);
                ++it2;
            }
        }
        c.append(InstructionType::TICK);
    }
    return std::move(c);
}
// Utility to partially merge two independent nodes, given a starting point
std::shared_ptr<CircuitNode> merge_nodes(std::shared_ptr<CircuitNode> nodea, std::shared_ptr<CircuitNode> nodeb, std::list<Instruction>::iterator indexa, std::list<Instruction>::iterator indexb)
{
    auto node0 = std::make_shared<CircuitNode>((nodea != nullptr ? nodea->name : "")+" + "+(nodeb != nullptr ? nodeb->name : ""));
    while (nodea != nullptr || nodeb != nullptr)
    {
        while (nodea != nullptr && indexa != nodea->circuit.instructions.end()) {
            if (indexa->type == InstructionType::TICK) {
                ++indexa;
                break;
            }
            node0->circuit.append(*indexa);
            ++indexa;
        }
        while (nodeb != nullptr && indexb != nodeb->circuit.instructions.end()) {
            if (indexb->type == InstructionType::TICK) {
                ++indexb;
                break;
            }
            node0->circuit.append(*indexb);
            ++indexb;
        }
        bool enda = nodea != nullptr && indexa == nodea->circuit.instructions.end();
        bool endb = nodeb != nullptr && indexb == nodeb->circuit.instructions.end();
        if (enda && endb) {
            if (!nodea->childs.empty() && !nodeb->childs.empty()) {
                for (int i=0; i<nodea->childs.size(); i++) {
                    for (int j=0; j<nodeb->childs.size(); j++) {
                        node0->childs.push_back(merge_nodes(nodea->childs[i], nodeb->childs[j], nodea->childs[i] == nullptr ? std::list<Instruction>::iterator() : nodea->childs[i]->circuit.instructions.begin(), nodeb->childs[j] == nullptr ? std::list<Instruction>::iterator() : nodeb->childs[j]->circuit.instructions.begin()));
                    }
                }
                node0->next_node_index = [nodea, nodeb](MeasurementResults &results) {
                    int i = nodea->next_node_index != nullptr ? nodea->next_node_index(results) : 0;
                    int j = nodeb->next_node_index != nullptr ? nodeb->next_node_index(results) : 0;
                    if (i < 0 || j < 0)
                        return -1;
                    return (int)(i*nodeb->childs.size()+j);
                };
            } else if (!nodea->childs.empty()) {
                node0->childs = nodea->childs;
                node0->next_node_index = nodea->next_node_index;
            } else if (!nodeb->childs.empty()) {
                node0->childs = nodeb->childs;
                node0->next_node_index = nodeb->next_node_index;
            }
            if (nodea->error_corrections && nodeb->error_corrections) {
                node0->error_corrections = [nodea, nodeb](MeasurementResults &results) {
                    auto errors = nodea->error_corrections(results);
                    auto errb = nodeb->error_corrections(results);
                    for (int err : errb.first) {
                        if (!errors.first.erase(err))
                            errors.first.insert(err);
                    }
                    for (int err : errb.second) {
                        if (!errors.second.erase(err))
                            errors.second.insert(err);
                    }
                    return errors;
                };
            } else if (nodea->error_corrections) {
                node0->error_corrections = nodea->error_corrections;
            } else {
                node0->error_corrections = nodeb->error_corrections;
            }
            break;
        } else if (enda && !nodea->childs.empty()) {
            if (nodeb == nullptr) {
                node0->childs = nodea->childs;
            } else {
                for (auto node : nodea->childs) {
                    node0->childs.push_back(merge_nodes(node, nodeb, node == nullptr ? std::list<Instruction>::iterator() : node->circuit.instructions.begin(), indexb));
                }
            }
            node0->next_node_index = nodea->next_node_index;
            node0->error_corrections = nodea->error_corrections;
            break;
        } else if (endb && !nodeb->childs.empty()) {
            if (nodea == nullptr) {
                node0->childs = nodeb->childs;
            } else {
                for (auto node : nodeb->childs) {
                    node0->childs.push_back(merge_nodes(nodea, node, indexa, node == nullptr ? std::list<Instruction>::iterator() : node->circuit.instructions.begin()));
                }
            }
            node0->next_node_index = nodeb->next_node_index;
            node0->error_corrections = nodeb->error_corrections;
            break;
        } else if (enda && nodeb == nullptr) {
            node0->circuit.append(InstructionType::TICK);
            node0->next_node_index = nodea->next_node_index;
            node0->error_corrections = nodea->error_corrections;
            break;
        } else if (endb && nodea == nullptr) {
            node0->circuit.append(InstructionType::TICK);
            node0->next_node_index = nodeb->next_node_index;
            node0->error_corrections = nodeb->error_corrections;
            break;
        }
        node0->circuit.append(InstructionType::TICK);
    }
    return node0;
}
// Combine two trees into a single one. Not working if the tree contains a loop
std::shared_ptr<CircuitNode> merge_nodes(std::shared_ptr<CircuitNode> nodea, std::shared_ptr<CircuitNode> nodeb)
{
    return merge_nodes(nodea, nodeb, nodea->circuit.instructions.begin(), nodeb->circuit.instructions.begin());
}
// Not directly called, called from the main apply_node_to_end to properly handle tree loops
void apply_node_to_end(std::shared_ptr<CircuitNode> node0, std::shared_ptr<CircuitNode> node, std::set<std::shared_ptr<CircuitNode>> &visited, std::shared_ptr<CircuitNode> ft_node)
{
    if (visited.find(node0) != visited.end())
        return;
    visited.insert(node0);
    if (node0->childs.empty()) {
        if (ft_node != nullptr)
            node0->childs.push_back(ft_node);
    } else {
        for (int i=0; i<node0->childs.size(); ++i) {
            if (i > 0)
                ft_node = node;
            if (node0->childs[i] == nullptr) {
                node0->childs[i] = ft_node;
            } else {
                apply_node_to_end(node0->childs[i], node, visited, ft_node);
            }
        }
    }
}
// Appends a tree to the end of another one. It is possible to append two different trees, one for the FT
// path (the one followed if no error happens) and another one for the rest of (non-FT) paths.
// This can be used, e.g., to switch to non-FT gadgets once an error has been detected.
void apply_node_to_end(std::shared_ptr<CircuitNode> node0, std::shared_ptr<CircuitNode> node, std::shared_ptr<CircuitNode> ft_node)
{
    std::set<std::shared_ptr<CircuitNode>> visited;
    apply_node_to_end(node0, node, visited, ft_node);
}
#include <sstream>
// Writes all possible circuit paths that a tree can contain
void print_nodes(std::shared_ptr<CircuitNode> node, std::set<std::shared_ptr<CircuitNode>> &&visited, std::string printed)
{
    if (printed.empty()) printed += "Branch 0 - ";
    printed += node->name+"\n";
    std::stringstream s;
    s<<node->circuit;
    printed += s.str();
    if (node->childs.empty()) {
        std::cout<<printed<<std::endl;
        return;
    }
    for (int i=0; i<node->childs.size(); i++) {
        if (node->childs[i]) {
            if (visited.find(node->childs[i]) == visited.end()) {
                std::set<std::shared_ptr<CircuitNode>> v = visited;
                v.insert(node->childs[i]);
                print_nodes(node->childs[i], std::move(v), printed+"Branch "+std::to_string(i)+" - ");
            } else {
                printed += "Branch "+std::to_string(i)+" - go back to "+node->childs[i]->name;
            }
        } else {
            //std::cout<<printed+std::to_string(i)<<" - END\n"<<std::endl;
        }
    }
}
std::vector<int> node_depth(std::shared_ptr<CircuitNode> node0)
{
    std::vector<int> depth0 = {0};
    for (int i=0; i<node0->childs.size(); i++) {
        if (node0->childs[i]) {
            auto depth = node_depth(node0->childs[i]);
            for (int j=0; j<depth.size(); j++) {
                if (j+1<depth0.size())
                    depth0[j+1] += depth[j];
                else
                    depth0.push_back(depth[j]);
            }
            depth0[0] += 1;
        }
    }
    return depth0;
}
std::set<CircuitNode*> visited;
int node_count(std::shared_ptr<CircuitNode> node0)
{
    visited.insert(node0.get());
    int count = 1;
    for (auto node : node0->childs) {
        if (node != nullptr)
            count += node_count(node);
    }
    return count;
}
void cnot_count(std::shared_ptr<CircuitNode> node, std::set<std::shared_ptr<CircuitNode>> &&visited, int current_count)
{
    for (auto &inst : node->circuit.instructions) {
        if (inst.type == InstructionType::CX)
            current_count += inst.targets.size()/2;
    }
    if (node->childs.empty()) {
        std::cout<<current_count<<std::endl;
        return;
    }
    for (int i=0; i<node->childs.size(); i++) {
        if (node->childs[i]) {
            if (visited.find(node->childs[i]) == visited.end()) {
                std::set<std::shared_ptr<CircuitNode>> v = visited;
                v.insert(node->childs[i]);
                cnot_count(node->childs[i], std::move(v), current_count);
            }
        }
    }
}