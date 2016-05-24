#' enaMTI --- Mixed Trophic Impacts Analysis
#' follows Ulanowicz and Puccia, 1990.
#' INPUT = network object
#' OUTPUT = list of trophic impact statistics
#' Borrett | June 2012, MKL | July 2013
#' ------------------------------------

enaMTI <- function(x,eigen.check=TRUE,zero.na=TRUE, balance.override=FALSE){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
                                        #Data checks
  if (any(is.na(x%v%'respiration'))){
    G <- FP <- Q <- M <- as.matrix(x, attrname = 'flow')
    G[is.na(G)==FALSE] <- FP[is.na(FP)==FALSE] <- Q[is.na(Q)==FALSE] <- M[is.na(M)==FALSE] <- NA
    out <- list('G'=G,'FP'=FP,'Q'=Q,'M'=M)
    warning('Model is missing respiration. Output is NA.')
  }else{
                                        #Check for balancing
    if (balance.override){}else{
      if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
      if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
    }
                                        #Unpack
    Flow <- as.matrix(x, attrname = 'flow')  #flows
    input <- x%v%'input' #inputs
    output <- x%v%'output'
    resp <- x%v%'respiration'
    stor <- x%v%'storage' #storage values
    I <- Flow*0;
    diag(I) <- 1 #create the identity matrix
    T. <- input + apply(Flow,2,sum)
                                        #
    G <- t(t(Flow)/T.)        # this is now oriented as in ulanowicz and puccia row to column
    FP <- Flow / (T.-resp)    # Authors exclude respiration to divide only by the NET production of the compartment.

                                        # check and replace NA values with 0 if zero.na
    if (zero.na){
      G[is.na(G)] <- 0
      FP[is.na(FP)] <- 0
    }
                                        #Make infinity values equal to zero
  G[is.infinite(G)] <- 0
  FP[is.infinite(FP)] <- 0
                                        # Set FP to zero when receiver compartment (j) is non-living
    FP[,which(x%v%'living'==FALSE)] <- 0
    Q <- G- t(FP)
    dom1Q <- abs(eigen(Q)$values[1])
    if(dom1Q <= 1 ){
      M <- ginv(I-Q)-I              # Total Impacts of i on j.
    } else {
      if(eigen.check==FALSE){
        M <- ginv(I-Q)-I              # Total Impacts of i on j.
      } else { M <- NA}
    }
    out <- list('G'=G,'FP'=FP,'Q'=Q,'M'=M)
  }
    return(out)
}
