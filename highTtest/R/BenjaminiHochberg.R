#------------------------------------------------------------------------------#
# Benjamini Hochberg method                                                    #
#------------------------------------------------------------------------------#
#                                                                              #
#  Inputs :                                                                    #
#                                                                              #
# alphas  : vector of levels at which the FDR is to be 'controlled'            #
#                                                                              #
# p_value : vector of p-values                                                 # 
#                                                                              #
#  Outputs :                                                                   #
#                                                                              #
#  matrix indicating rejected hypotheses for each level                        #
#------------------------------------------------------------------------------#
BenjaminiHochberg <- function(alphas, 
                              p_value){

  adjusted.p_value <- p.adjust(p_value, method="BH")

  indicator.bh <- sapply(X = alphas, 
                         FUN = function(x, p){ p <= x }, 
                         p = adjusted.p_value)

  if(is.matrix(indicator.bh)){
    colnames(indicator.bh) <- round(alphas,3)
  } else {
    indicator.bh <- matrix(indicator.bh,nrow=1)
    colnames(indicator.bh) <- round(alphas,3)
  }

  return(indicator.bh)

}
