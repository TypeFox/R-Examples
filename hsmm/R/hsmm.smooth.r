hsmm.smooth <- function(x, 
                       od, 
                       rd, 
                       pi.par,
                       tpm.par,
                       od.par,
                       rd.par,
                       M = NA)
{
  
  # formatting, renaming, additional variables
  inputData <- as.vector(x)
  tau       <- get.tau(inputData)
  error     <- as.integer(0)

  # smoohting results are independent of censoring
  censoring <- as.integer(0)
   
  # determine number of states
  J <- length(pi.par)

  # write all initial values into one list
  Para     <- list()
  Para$pi  <- pi.par 
  Para$tpm <- t(tpm.par) 
  Para$rd  <- rd.par 
  Para$od  <- od.par 

  # selection of the maximum runlength
  if (is.na(M)){
    if (rd == "nonp"){
      M <- as.integer(dim(rd.par$np)[1])
      } else {
      M <- as.integer(min(tau, 1000))
      }
    } # endif isna(M)
  
  # variables calculated by Forward-Backward alg.
  F    <- as.double(rep(0, times = J * tau))
  L    <- as.double(rep(0, times = J * tau))
  G    <- as.double(rep(0, times = J * tau))
  L1   <- as.double(rep(0, times = J * tau))
  N    <- as.double(rep(0, times = tau))
  Norm <- as.double(rep(0, times = J * tau))
  eta  <- as.double(rep(0, times = J * M))
  xi   <- as.double(rep(0, times = J * M))

  # initialize variables for parameters estimated by EM
  EM.Para <- Para
    
  # Store variables for calling FB
  FB.p.tpm      <- EM.Para$tpm
  dim(FB.p.tpm) <- c(J * J)
  FB.pi.ini     <- EM.Para$pi
  FB.d          <- get.d(rd, J, M, param = EM.Para$rd)
  FB.pdf        <- get.pdf(inputData, od, J, M, param = EM.Para$od)   

  # Call Forward-Backward Alg.
  FB.result <- FB(censoring, tau, J, M, FB.d, FB.p.tpm, FB.pi.ini, FB.pdf, F, L, G, L1, N, Norm, eta, xi, error)
  error     <- FB.result[[17]]
    
  if (error == 0){   
    # Save results of FB
    F    <- FB.result[[9]]
    L    <- FB.result[[10]]
    G    <- FB.result[[11]]
    L1   <- FB.result[[12]]
    N    <- FB.result[[13]]
    Norm <- FB.result[[14]]
    eta  <- FB.result[[15]]
    xi   <- FB.result[[16]]
    
    # Change results of FB to matrices
    dim(F)    <- c(tau, J)
    dim(L)    <- c(tau, J)
    dim(G)    <- c(tau, J)
    dim(L1)   <- c(tau, J)
    dim(Norm) <- c(tau, J)
    dim(eta)  <- c(M, J)
    dim(xi)   <- c(M, J)
    
    # Calculate llh
    llh <- sum(log(N[1:tau]))
    
    } # Has an error occured?  
  
  # return results
  out <- list(call        = match.call(),
              smooth.prob = t(L),
              path        = apply(t(L), 2, which.max))
  return(out)
  }
  