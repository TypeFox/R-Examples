#------------------------------------------------------------------------------#
# Storey Tibshirani method                                                     #
# This code is based on that provided by John Storey obtained from             #
# http://genomics.princeton.edu/storeylab/qvalue/results.R                     #
#------------------------------------------------------------------------------#
#                                                                              #
#  Inputs :                                                                    #
#                                                                              #
#  alphas : vector of levels at which the FDR is to be 'controlled'            #
#                                                                              #
# p_value : vector of p-values                                                 # 
#                                                                              #
#  Outputs :                                                                   #
#                                                                              #
#  logical matrix indicating rejected hypotheses for each level                #
#------------------------------------------------------------------------------#
StoreyTibshirani <- function(alphas, p_value){

  m <- length(p_value)

  lam <- seq(0,0.95,0.01)

  pi0_temp <- sapply(X = lam,
                     FUN = function(x,p){ mean(p>x)/(1-x)},
                     p = p_value)

  pi0 <- smooth.spline(x = lam, 
                       y = pi0_temp, 
                       w = (1-lam), 
                       df = 3)$y[length(lam)]

  qvalue <- pi0*m*p_value/rank(p_value)

  qvalue[m] <- min(qvalue[m],1)

  i <- m-1
  while(i >= 1) {
    qvalue[i] <- min(qvalue[i],qvalue[i+1],1)
    i <- i - 1
  }

  indicator.st <- sapply(X = alphas, 
                         FUN = function(x,q){q <= x}, 
                         q = qvalue)

  if(is.matrix(indicator.st)){
    colnames(indicator.st) <- round(alphas,3)
  } else {
    indicator.st <- matrix(indicator.st,nrow=1)
    colnames(indicator.st) <- round(alphas,3)
  }

  return(indicator.st)
}
