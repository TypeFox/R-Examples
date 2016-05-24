# Define Interface to C++
Viterbi <- function(tau, J, M, VT.d, VT.p.tpm, VT.pi.ini, VT.pdf, VT.hiddenStates){
  .C("Viterbi", tau, J, M, VT.d, VT.p.tpm, VT.pi.ini, VT.pdf, VT.hiddenStates, PACKAGE="hsmm")
  }    
