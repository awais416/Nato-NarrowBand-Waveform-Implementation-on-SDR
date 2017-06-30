## Summary
To design an optimum transceiver using Continuous Phase Modulation (CPM).  The bandwidth given in specification was 25kH and data rate to achieve was up to 96Kbps on synthetic data. Transmitter Involved quantization, encryption, and modulation. Complete Matlab simulation of a transceiver. Programmed transmitter to C Language for DSP and receiver in Verilog for FPGA.  An optimum receiver designed using Viterbi algorithm &amp; MLSE. Involved decryption and demod. 

**Matlab Simulation**
Main code.m contains matlab simulation which is a complete transceiver iteself. It was made to analyze how data rate is effected by applying different algorithms of modulation and encryption

**main.c**
This is a complete transmitter which was programmed in C language by converting the transmitter part from Main code.m file. This code was further burned in DSP kit of SDR for transmitting

**Narrowband-verilog_receiver.rar**
The receiver was converted into verilog language from the Main code.m matlab file. It was further burned in FPGA kit of SDR for receiving data.
