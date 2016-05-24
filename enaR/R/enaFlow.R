#' enaFlow --- flow analysis
#' INPUT = network object
#' OUTPUT = list of flow statistics
#'
#' M. Lau | July 2011
#' ---------------------------------------------------

enaFlow <- function(x,zero.na=TRUE,balance.override=FALSE){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }

                                        # unpack model
  Flow <- t(as.matrix(x,attrname = 'flow')) #flows
  input <- x%v%'input' #inputs
  stor <- x%v%'storage' #storage values

  n <- nrow(Flow)      # number of nodes
  I <- diag(1,nrow(Flow),ncol(Flow))          # create identity matrix
  T. <- apply(Flow,1,sum) + input;   # input throughflow (assuming steady state)

                                        #compute the intercompartmental flows
  GP <- Flow / T.  #Input perspective
  G <- t(t(Flow) / T.)  #Output perspective

                                        #check and replace NA values with 0 if zero.na
  if (zero.na){
    GP[is.na(GP)] <- 0
    G[is.na(G)] <- 0
    GP[is.infinite(GP)] <- 0
    G[is.infinite(G)] <- 0
  }

                                        #compute the integral flows
  NP <- ginv((I - GP))
  rownames(NP) <- colnames(NP) <- colnames(GP)
  N <- ginv((I - G))
  rownames(N) <- colnames(N) <- colnames(G)

  ## the ginv function creates noticible numeric error.  I am removing some of it here by rounding
  tol <- 10
  N <- round(N,tol)
  NP <- round(NP,tol)


  ## Network Statistics
  TST <- sum(T.)  # total system throughflow
  TSTp <- sum(Flow) + sum(x%v%'input') + sum(x%v%'output') # total system throughput

  Boundary <- sum(input)
  APL <- TST/Boundary  # Average Path Lenght (Finn 1976; aka network
                      # aggradation, multiplier effect)

                                        # Finn Cycling Index
  p <- as.matrix(rep(1,n),nrow=n)
  dN <- diag(N)
  TSTc <- sum((dN-p)/dN *T.)
  FCI <- TSTc/TST

                                        # non-locality (realized)
  direct <- sum(G %*% input)
  indirect <- sum((N - I - G) %*% input)
  ID.F <- indirect/direct
  BFI <- Boundary/TST
  DFI <- sum(G %*% input) / TST
  IFI <- indirect/TST
                                        # non-locality (idealized)
  ID.F.O <- sum(N-I-G)/sum(G)
  ID.F.I <- sum(NP-I-GP)/sum(GP)
                                        # HMG
  HMG.O <- ( sd(as.vector(G)) / mean(G) ) / ( sd(as.vector(N)) / mean(N) )
  HMG.I <- ( sd(as.vector(GP)) / mean(GP) ) / ( sd(as.vector(NP)) / mean(NP) )
                                        # Amplification
  AMP.O <- length(which( (N - diag(diag(N))) > 1))
  AMP.I <- length(which( (NP - diag(diag(NP))) > 1))
                                        # MODE ANALYSIS
                                        # This is built from Fath's original MODE program
  mode0.F <- Boundary                     # boundary flow
  mode1.F <- sum(ginv(diag(diag(N))) %*% N %*% diag(input) - diag(as.vector(I%*%input))) # internal first passage flow
  mode2.F <- sum((diag(diag(N))-I) %*% ginv(diag(diag(N))) %*% N %*% diag(input))  # cycled flow
  mode3.F <- mode1.F                      # dissipative equivalent to mode 1
  mode4.F <- mode0.F                    # dissipative equivalent to mode 0
                                        #re-orientation
  orient <- get.orient()
  if (orient == 'rc'){
    G <- t(G)
    GP <- t(GP)
    N <- t(N)
    NP <- t(NP)
  }else{}
                                        #network statistics
  ns <- cbind(Boundary,TST,TSTp,APL,FCI,
              BFI,DFI,IFI,
              ID.F,ID.F.I,ID.F.O,
              HMG.I,HMG.O,
              AMP.I,AMP.O,
              mode0.F,mode1.F,mode2.F,mode3.F,mode4.F)

                                        #output
  return(list('T'=T.,'G'=G,'GP'=GP,'N'=N,'NP'=NP,'ns'=ns))
}

