#' enaAscendency --- calculates the ascendency statistics
#' of Ulanowicz
#' INPUT = network object
#' OUTPUT = matrix of ascendency statistics
#' 
#' D. Hines | December 2011
#' ---------------------------------------------------

enaAscendency <- function(x='network object'){

  if (any(is.na(x%v%'export'))|any(is.na(x%v%'respiration'))){
    warning('Model is missing either export or respiration. Calculations may not be correct.')
  }

#'####### set initial conditions for calculations #########
  T.ulan <- as.extended(x)
  N <- ncol(T.ulan) # set up N
  r.td <- c.ld <- t.ulan <- ami <- mat.or.vec(N,N) # initialize ascendency matrix
  oh <- mat.or.vec(N,N) # initialize overhead matrix
  cap <- mat.or.vec(N,N) # initialize capacity matrix
                                        #calculate total system throughPUT
  TSTp <- sum(T.ulan)
  
#'################### calculate AMI #######################
  
                                        # loop through T.ulan to calculate AMI
  for (i in 1:N){
    for (j in 1:N){
      if (T.ulan[i,j] == 0){
        ami[i,j] <- 0
      }else{
        ami[i,j] <- T.ulan[i,j]/TSTp * log2((T.ulan[i,j]*TSTp)/(sum(T.ulan[,j])*sum(T.ulan[i,])))
      }
    }
  }
  
  AMI <- sum(ami)
  
#'################ calculate ascendency ###################
  
  ASC <- TSTp*AMI
  
#'################ calculate overhead  ####################
  
                                        # loop through T.ulan to calculate overhead
  for (i in 1:N){
    for (j in 1:N){
      if (T.ulan[i,j] == 0){
        oh[i,j] <- 0
      }else{
        oh[i,j] <- T.ulan[i,j] * log2((T.ulan[i,j]^2)/(sum(T.ulan[,j])*sum(T.ulan[i,])))
      }
    }
  }
  
  OH <- -sum(oh)
  
#'############## calculate capacity (long) ###############
  
                                        # loop through T.ulan to calculate capacity
  for (i in 1:N){
    for (j in 1:N){
      if (T.ulan[i,j] == 0){
        cap[i,j] <- 0
      }else{
        cap[i,j] <- T.ulan[i,j] * log2((T.ulan[i,j])/(TSTp))
      }
    }
  }
  
  CAP <- -sum(cap)
  
#'############# calculate capacity (short) ################
  
  CAP2 <- ASC+OH
  
  
#'################### calculate ratios ####################
  
                                        # ratio for ascendency/capacity
  ASC.CAP <- ASC/CAP
  
                                        #ratio for overhead/capacity
  OH.CAP <- OH/CAP
  
                                        #confirm ratios sum to 1
  ASC.OH.RSUM <- ASC.CAP + OH.CAP

  robustness = -1 * ASC.CAP * log(ASC.CAP)  # robustness from Ulanowicz 2009; Fath 2014
  
  ################# Calculating Effective Link Density and Trophic Depth ########
  ## Calculate t.ulan 't'
  
  for (i in 1:N) {
        for (j in 1:N) {
            if (T.ulan[i, j] == 0) {
                t.ulan[i, j] <- 0
            }
            else {
                t.ulan[i, j] <- T.ulan[i,j]/TSTp
            }
        }
    }
    
    ## Effective Link Density (c)
    for (i in 1:N) {
        for (j in 1:N) {
            if (t.ulan[i, j] == 0) {
                c.ld[i, j] <- 1
            }
            else {
                c.ld[i, j] <- (sqrt(sum(t.ulan[i,])*sum(t.ulan[,j]))/t.ulan[i,j])^(t.ulan[i,j])
            }
        }
    }
    C.LD <- prod(c.ld)
    
    ## Trophic Depth (r)
    for (i in 1:N) {
        for (j in 1:N) {
            if (t.ulan[i, j] == 0) {
                r.td[i, j] <- 1
            }
            else {
                r.td[i, j] <- (t.ulan[i,j]/(sum(t.ulan[i,])*sum(t.ulan[,j])))^(t.ulan[i,j])
            }
        }
    }
    R.TD <- prod(r.td)
    
    ELD <- C.LD
    TD <- R.TD
    ##############################################################################
    
    ns <- cbind(AMI,ASC,OH,CAP,ASC.CAP,OH.CAP, robustness, ELD, TD)
                                        #
  return(ns)
  
}
