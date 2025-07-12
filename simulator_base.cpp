#include "simulator_base.h"
// Introduces a X error with probability p
void BaseFrameSimulator::x_error(int qubit, double p)
{
    std::geometric_distribution<size_t> dist(p == 1 ? 0.5 : p);
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (p == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        std::cerr<<"X"<<qubit<<" at tick "<<current_tick<<std::endl;
#endif
        flip_error(shot, qubit, ERROR_X);
    }
}
// Introduces a Y error with probability p
void BaseFrameSimulator::y_error(int qubit, double p)
{
    std::geometric_distribution<size_t> dist(p == 1 ? 0.5 : p);
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (p == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        std::cerr<<"Y"<<qubit<<" at tick "<<current_tick<<std::endl;
#endif
        flip_error(shot, qubit, ERROR_X|ERROR_Z);
    }
}
// Introduces a Z error with probability p
void BaseFrameSimulator::z_error(int qubit, double p)
{
    std::geometric_distribution<size_t> dist(p == 1 ? 0.5 : p);
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (p == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        std::cerr<<"Z"<<qubit<<" at tick "<<current_tick<<std::endl;
#endif
        flip_error(shot, qubit, ERROR_Z);
    }
}
void BaseFrameSimulator::depolarize(std::vector<int> &qubits, double p)
{
    std::geometric_distribution<size_t> dist(p == 1 ? 0.5 : p);
    std::uniform_int_distribution<> depol(1, (4<<(2*qubits.size()-2))-1);
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (p == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
        int type = depol(rng);
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        /*std::string names[] = {"","X","Z","Y"};
        std::cout<<names[type]<<qubit<<" TICK "<<current_tick<<std::endl;*/
#endif
        for (int i=0; i<qubits.size(); i++) {
            flip_error(shot, qubits[i], (type>>(2*i))&3);
        }
    }
}
// Introduces a depolarizing error with probability p
// If error happens, randomly choose between X, Y or Z errors
void BaseFrameSimulator::depolarize1(int qubit, double p)
{
    std::geometric_distribution<size_t> dist(p == 1 ? 0.5 : p);
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (p == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
        int type = depolarizer1(rng);
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        std::string names[] = {"","X","Z","Y"};
        std::cerr<<names[type]<<qubit<<" at tick "<<current_tick<<std::endl;
#endif
        flip_error(shot, qubit, type);
    }
}
// Introduces a two-qubit depolarizing error with probability p
// If error happens, randomly choose between {I,X,Y,Z}^2-{IxI} errors
void BaseFrameSimulator::depolarize2(int control, int target, double p)
{
    std::geometric_distribution<size_t> dist(p == 1 ? 0.5 : p);
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (p == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
        int type = depolarizer2(rng);
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        std::string names[] = {"","X","Z","Y"};
        std::cerr<<names[type&3]<<control<<"*"<<names[type>>2]<<target<<" at tick "<<current_tick<<std::endl;
#endif
        flip_error(shot, control, type & 3);
        flip_error(shot, target, type>>2);
    }
}
void BaseFrameSimulator::pauli1(int qubit, std::vector<double> &p)
{
    double ptot = 0;
    for (int i=0; i<3; i++) {
        ptot += p[i];
    }
    std::geometric_distribution<size_t> dist(ptot == 1 ? 0.5 : ptot);
    std::discrete_distribution<int> dist2(p.begin(), p.end());
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (ptot == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
        int type = dist2(rng)+1;
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        std::string names[] = {"I","X","Z","Y"};
        std::cerr<<names[type]<<qubit<<" at tick "<<current_tick<<std::endl;
#endif
        flip_error(shot, qubit, type);
    }
}
void BaseFrameSimulator::pauli2(int control, int target, std::vector<double> &p)
{
    double ptot = 0;
    for (int i=0; i<3; i++) {
        ptot += p[i];
    }
    std::geometric_distribution<size_t> dist(ptot == 1 ? 0.5 : ptot);
    std::discrete_distribution<int> dist2(p.begin(), p.end());
    size_t next_candidate = 0;
    for (;;) {
        size_t shot = next_candidate + (ptot == 1 ? 0 : dist(rng));
        if (shot >= num_shots)
            break;
        next_candidate = shot+1;
        int type = dist2(rng)+1;
#ifdef CHECK_FT
        if (error)
            continue;
        error = true;
        std::string names[] = {"","X","Z","Y"};
        std::cerr<<names[type&3]<<control<<"*"<<names[type>>2]<<target<<" at tick "<<current_tick<<std::endl;
#endif
        flip_error(shot, control, type & 3);
        flip_error(shot, target, type>>2);
    }
}
// Runs a circuit instruction
void BaseFrameSimulator::run(Instruction &instruction)
{
    //std::cout<<instruction<<" TICK "<<current_tick<<std::endl;
    switch (instruction.type) {
        case InstructionType::CX:
            for (int i=0; i<instruction.targets.size(); i+=2)
                cx(instruction.targets[i], instruction.targets[i+1]);
            break;
        case InstructionType::CZ:
            for (int i=0; i<instruction.targets.size(); i+=2)
                cz(instruction.targets[i], instruction.targets[i+1]);
            break;
        case InstructionType::SXX:
        case InstructionType::SXXDG:
            for (int j=1; j<instruction.targets.size(); j++) {
                for (int i=0; i<j; i++) {
                    sxx(instruction.targets[i], instruction.targets[j]);
                }
            }
            break;
        case InstructionType::SZZ:
        case InstructionType::SZZDG:
            for (int j=1; j<instruction.targets.size(); j++) {
                for (int i=0; i<j; i++) {
                    szz(instruction.targets[i], instruction.targets[j]);
                }
            }
            break;
        case InstructionType::MX:
            for (int i=0; i<instruction.targets.size(); i++) {
                mx(instruction.targets[i], *instruction.measurement_tag);
            }
            break;
        case InstructionType::MZ:
            for (int i=0; i<instruction.targets.size(); i++) {
                mz(instruction.targets[i], *instruction.measurement_tag);
            }
            break;
        case InstructionType::RX:
            for (int i=0; i<instruction.targets.size(); i++)
                rx(instruction.targets[i]);
            break;
        case InstructionType::RZ:
            for (int i=0; i<instruction.targets.size(); i++)
                rz(instruction.targets[i]);
            break;
        case InstructionType::H:
        case InstructionType::SY:
        case InstructionType::SYDG:
            for (int i=0; i<instruction.targets.size(); i++)
                h(instruction.targets[i]);
            break;
        case InstructionType::S:
        case InstructionType::SDG:
            for (int i=0; i<instruction.targets.size(); i++)
                s(instruction.targets[i]);
            break;
        case InstructionType::SX:
        case InstructionType::SXDG:
            for (int i=0; i<instruction.targets.size(); i++)
                sx(instruction.targets[i]);
            break;
        case InstructionType::X_ERROR: 
            for (int i=0; i<instruction.targets.size(); i++)
                x_error(instruction.targets[i], instruction.p[0]);
            break;
        case InstructionType::Y_ERROR: 
            for (int i=0; i<instruction.targets.size(); i++) {
                y_error(instruction.targets[i], instruction.p[0]);
            }
            break;
        case InstructionType::Z_ERROR: 
            for (int i=0; i<instruction.targets.size(); i++)
                z_error(instruction.targets[i], instruction.p[0]);
            break;
        case InstructionType::DEPOLARIZE: 
            depolarize(instruction.targets, instruction.p[0]);
            break;
        case InstructionType::DEPOLARIZE1: 
            for (int i=0; i<instruction.targets.size(); i++)
                depolarize1(instruction.targets[i], instruction.p[0]);
            break;
        case InstructionType::DEPOLARIZE2: 
            for (int i=0; i<instruction.targets.size(); i+=2)
                depolarize2(instruction.targets[i], instruction.targets[i+1], instruction.p[0]);
            break;
        case InstructionType::PAULI1:
            for (int i=0; i<instruction.targets.size(); i++)
                pauli1(instruction.targets[i], instruction.p);
            break;
        case InstructionType::PAULI2: 
            for (int i=0; i<instruction.targets.size(); i+=2)
                pauli2(instruction.targets[i], instruction.targets[i+1], instruction.p);
            break;
        case InstructionType::TICK:
            ++current_tick;
            break;
    }
}
// Runs a deterministic circuit
void BaseFrameSimulator::run(Circuit &circuit)
{
    for (auto &i : circuit.instructions) {
        run(i);
    }
}