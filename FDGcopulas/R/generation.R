
## Generation 
# - n : sample size
# - FDGcopula : the FDGcopula object
# - sizeSubSample : size of the sample the maxima are taken over in
# the generation of extreme-value copulas

rFDG <- function(n, FDGcopula, sizeSubSample=10000){
  famId <- switch(FDGcopula@family,
                  "frechet" = 1,
                  "cuadrasauge" = 2,
                  "sinus" = 3,
                  "exponential" = 4)
  if(FDGcopula@extremevalue){
    randGen_ev_CPP_2(n, famId, FDGcopula@parameters, sizeSubSample)
  }else{
    randGenCPP_2(n, famId, FDGcopula@parameters)
  }
}
