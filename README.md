# Decoding-HMM-for-Molecular-Sequences
This codes aim to find the path through the HMM that maximizes the loglikelihood
of our observed molecular sequence x = x1, x2, . . . , xL, that is,
$y^∗ = arg max (log P[X = x, Y = y|θ])$

This can beachieved using the Viterbi algorithm in HMM.py code.
In The Viterbi algorithm, the best path y∗ can be recovered by backtracking through the V “matrix”,
which is typically implemented by storing additional information about the best paths for submodels
when filling in the V matrix.
