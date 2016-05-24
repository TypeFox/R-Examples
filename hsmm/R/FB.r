# Define Interface to C++
FB <- function(censoring, tau, J, M, FB.d, FB.p.tpm, FB.pi.ini, FB.pdf, 
               F, L, G, L1, N, Norm, eta, xi, error){
  .C("FB", censoring, tau, J, M, FB.d, FB.p.tpm, FB.pi.ini, FB.pdf, 
     F, L, G, L1, N, Norm, eta, xi, error, PACKAGE="hsmm");
  }
 