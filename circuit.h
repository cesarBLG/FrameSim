#pragma once
#include <vector>
#include <memory>
#include <list>
#include <set>
#include <functional>
#include <map>
#include <ostream>
#include <optional>
enum struct InstructionType
{
    I,
    X,
    Y,
    Z,
    H,
    S,
    SDG,
    SX,
    SXDG,
    SY,
    SYDG,
    CX,
    CY,
    CZ,
    SXX,
    SXXDG,
    SZZ,
    SZZDG,
    MX,
    MY,
    MZ,
    RX,
    RY,
    RZ,
    DEPOLARIZE,
    DEPOLARIZE1,
    DEPOLARIZE2,
    X_ERROR,
    Y_ERROR,
    Z_ERROR,
    PAULI1,
    PAULI2,
    DELAY,
    TICK,
};
// Measurement tag: allows identifying a measurement by its name and round
// A measurement tag must be provided for all measurements in a circuit,
// to properly identify it
struct MeasurementTag
{
    int current_round;
    std::string name;
    bool operator<(const MeasurementTag &o) const
    {
        if (current_round != o.current_round) return current_round < o.current_round;
        return name < o.name;
    }
};
// Contains measurement results from a specific shot, indicating whether it
// has been flipped with respect to the noiseless expected outcome due to an error
struct MeasurementResults
{
    virtual bool is_flipped(int qubit, const MeasurementTag &tag) = 0;
    virtual bool reset_flipped(int qubit, const MeasurementTag &tag) = 0;
    virtual void flip(int qubit, const MeasurementTag &tag) = 0;
};
// Instruction: corresponds to one or several gates of the same type applied to a set of qubits
struct Instruction
{
    InstructionType type;
    std::vector<int> targets;
    std::vector<double> p;
    std::optional<MeasurementTag> measurement_tag;
    std::optional<std::string> label;
    Instruction() = default;
    Instruction(InstructionType type, std::vector<int> targets) : type(type), targets(targets)
    {

    }
    Instruction(InstructionType type, std::vector<int> targets, double p) : type(type), targets(targets), p({p})
    {

    }
    Instruction(InstructionType type, std::vector<int> targets, std::vector<double> p, std::optional<MeasurementTag> tag=std::nullopt, std::string label="") : type(type), targets(targets), p(p), measurement_tag(tag), label(label)
    {

    }
};
// Represents a set of different quantum gates
class Circuit
{
    public:
    Circuit() : num_qubits(0) {}
    std::list<Instruction> instructions;
    int num_qubits;
    Circuit &operator += (const Circuit &o)
    {
        if (o.num_qubits > num_qubits)
            num_qubits = o.num_qubits;
        instructions.insert(instructions.end(), o.instructions.begin(), o.instructions.end());
        return *this;
    }
    Circuit operator+(const Circuit &o) const
    {
        Circuit c = *this;
        c += o;
        return std::move(c);
    }
    Circuit &append(InstructionType type)
    {
        return append(Instruction(type, {}));
    }
    Circuit &append(InstructionType type, int target1)
    {
        return append(Instruction(type, {target1}));
    }
    Circuit &append(InstructionType type, int target1, int target2)
    {
        return append(Instruction(type, {target1, target2}));
    }
    Circuit &append(Instruction inst)
    {
        for (auto t : inst.targets) {
            if (t+1 > num_qubits)
                num_qubits = t+1;
        }
        instructions.push_back(inst);
        return *this;
    }

};
// Represents a circuit which is followed by different circuits depending on previous measurement outcomes
// Used for in-sequence logic
// Each circuit node can have several nodes as successors
// Which node is selected from all childs is determined by the next_node_index() function,
// taking into account measurement history, which must be provided by the user
// Error correction can be applied by providing the error_correction() function,
// which should return a list of qubits to be flipped depending on the
// measurement outcomes.
class CircuitNode : public std::enable_shared_from_this<CircuitNode>
{
    public:
    std::string name;
    std::vector<std::shared_ptr<CircuitNode>> childs;
    Circuit circuit;
    // Determines next node depending on measurement outcomes
    // Parameter: list of flipped measurements, ordered by qubit
    // Must return index of next circuit, -1 if shot discarded by postselection
    std::function<int(MeasurementResults&)> next_node_index;
    // Applies error correction depending on measurement outcomes
    // Parameter: list of flipped measurements, ordered by qubit
    // Must return a qubits for which X errors (pair.first) and Z errors (pair.second) have to be applied to recover
    std::function<std::pair<std::set<int>,std::set<int>>(MeasurementResults&)> error_corrections;
    CircuitNode(std::string name) : name(name) {}
    std::shared_ptr<CircuitNode> deep_copy()
    {
        auto node = std::make_shared<CircuitNode>(name);
        node->circuit = circuit;
        node->next_node_index = next_node_index;
        node->error_corrections = error_corrections;
        for (auto &child : childs) {
            if (child)
                node->childs.push_back(child->deep_copy());
            else
                node->childs.push_back(nullptr);
        }
        return node;
    }
};
std::ostream &operator<<(std::ostream &os, const Instruction &instr);
std::ostream &operator<<(std::ostream &os, const Circuit &m);

Circuit merge_circuits(Circuit c1, Circuit c2);
std::shared_ptr<CircuitNode> merge_nodes(std::shared_ptr<CircuitNode> nodea, std::shared_ptr<CircuitNode> nodeb);
void apply_node_to_end(std::shared_ptr<CircuitNode> node0, std::shared_ptr<CircuitNode> node, std::shared_ptr<CircuitNode> ft_node);
void print_nodes(std::shared_ptr<CircuitNode> node, std::set<std::shared_ptr<CircuitNode>> &&visited = std::set<std::shared_ptr<CircuitNode>>(), std::string printed="");
std::vector<int> node_depth(std::shared_ptr<CircuitNode> node0);
int node_count(std::shared_ptr<CircuitNode> node0);
void cnot_count(std::shared_ptr<CircuitNode> node, std::set<std::shared_ptr<CircuitNode>> &&visited = std::set<std::shared_ptr<CircuitNode>>(), int current_count=0);