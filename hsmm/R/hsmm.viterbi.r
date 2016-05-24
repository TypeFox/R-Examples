hsmm.viterbi <- function(x, 
                         od, 
                         rd, 
                         pi.par,
                         tpm.par,
                         od.par,
                         rd.par,
                         M         = NA){
                    
  od.t <- c("binom", "norm", "pois", "t")
  rd.t <- c("nonp", "geom", "nbinom", "log", "pois")
  
  inputData <- as.vector(x)  
  tau       <- get.tau(inputData)
  # determine number of states
  J <- length(pi.par)

  # write all initial values into one list
  Para     <- list()
  Para$pi  <- pi.par 
  Para$tpm <- t(tpm.par) 
  Para$od  <- od.par 
  Para$rd  <- rd.par 

  # selection of the maximum runlength
  if (is.na(M)){
    if (rd == "nonp"){
      M <- as.integer(dim(rd.par$np)[1])
      } else {
      M <- as.integer(min(tau, 1000))
      }
    } # endif isna(M)
  
  # Store variables for calling FB
  VT.p.tpm      <- Para$tpm
  dim(VT.p.tpm) <- c(J * J)
  VT.pi.ini     <- Para$pi
  VT.d   <- get.d(rd, J, M, param = Para$rd)
  VT.pdf <- get.pdf(inputData, od, J, M, param = Para$od)
  VT.hiddenStates <- rep(as.integer(0), times = tau)
  
  # Call Viterbi
  VT.results <- Viterbi(tau, J, M, VT.d, VT.p.tpm, VT.pi.ini, VT.pdf, VT.hiddenStates)

  # Save results of Viterbi
  VT.hiddenStates <- VT.results[[8]] + 1

  out <- list(call = match.call(),
              path = VT.hiddenStates)
  return(out)
}
  
