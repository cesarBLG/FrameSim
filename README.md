# Sparse Frame Simulator
Sparse frame simulator, inspired by Stim, with support for in-sequence logic / dynamic circuits

By using a sparse representation of data, shots where no fault has yet ocurred have zero computational cost
Multiple shots are executed in parallel for more efficient noise sampling.

The basic objects in the simulator are:
- Instruction: indicates the clifford gate to be executed, the target qubits and gate arguments.
- Circuit: contains a list of instructions. Gates can be separated by TICK instructions.
- CircuitNode: contains a circuit, as well as a list with the possible nodes that go after the current one,
  defining a circuit tree. After a node is executed, the simulator will decide the next node to execute
  for every shot, depending on measurement outcomes. Error correction can be applied after each node.
- NoiseModel: user-defined class with a noisy_circuit() function that applies noise to a circuit or tree
- MeasurementTag: for every measurement in the circuit, it contains a string and integer identifying it
- MeasurementResults: records whether a measurement has been flipped due to noise. This is used for error
  correction and next-node determination, by providing custom functions to a CircuitNode

A dense frame simulator with support for dynamic circuits is also provided, to allow circuit validation.