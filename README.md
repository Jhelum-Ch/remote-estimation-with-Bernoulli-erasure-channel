# remote-estimation-with-Bernoulli-erasure-channel
This repository contains codes for remote estimation with two agents connected via a Bernoulli erasure channel. We have considered two models for the source, Model A and Model B.

Model A: A symmetric scalar Birth-Death Markov chain with countable state space; the birth rate is denoted by 'p'. From any state 'x', the transition probabilities to states 'x+1', 'x-1' and 'x' are 'p', 'p' and '1-2p' respectively. The code for this model is written in MATLAB. 

These codes reproduce the results presented in the paper:

J. Chakravorty and A. Mahajan, "Remote-state estimation with packet drop", in proceedings of the IFAC Workshop on Distributed Estimation and Control in Networked Systems, Tokyo, Japan, 2016.


Model B: First-order scalar autoregressive process with zer-mean Gaussian innovations. We have used simulation based method (we call it Renewal Monte Carlo) to solve the problem. The codes for this model are written in Julia. A documentation for the Renewal Monte Carlo is provided in a separate PDF file for ease of interpretation. 

These codes reproduce the results presented in the paper:

J. Chakravorty, J. Subramanian, A. Mahajan, "Stochastic approximation based methods for computing the optimal thresholds in remote-state estimation with packet drops," in proceedings of the American Control Conference, Seattle, WA, 2017.
