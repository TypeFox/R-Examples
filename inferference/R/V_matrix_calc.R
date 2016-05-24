#-----------------------------------------------------------------------------#
# V Matrix
# 
# Computes the V matrix necessary for variance estimates. Used in 
# \code{\link{ipw_effect_calc}}. See web appendix of Perez et al. 2014 for more details.
# 
# @param scores the output of \code{\link{score_matrix}}
# @param point_estimates output of \code{\link{ipw_point_estimates}}
# @param allocation1 See details in \code{\link{ipw_effect_calc}}.
# @param trt.lvl1 See details in \code{\link{ipw_effect_calc}}.
# @param allocation2 See details in \code{\link{ipw_effect_calc}}.
# @param trt.lvl2 See details in \code{\link{ipw_effect_calc}}.
# @param effect_type See details in \code{\link{ipw_effect_calc}}.
# @param marginal See details in \code{\link{ipw_effect_calc}}.
# @return V matrix
# @export
# 
#-----------------------------------------------------------------------------#
V_matrix <- function(scores, 
                     point_estimates, 
                     allocation1, 
                     trt.lvl1, 
                     allocation2 = NA, 
                     trt.lvl2    = NA, 
                     effect_type, 
                     marginal){
  ## Necessary bits ##
  N  <- dim(scores)[1]
  p  <- dim(scores)[2]
  a1 <- allocation1
  a2 <- allocation2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  ## Grab the last element of the psi(O, theta) vector: psi_a, alpha ##
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  hold_oal <- point_estimates[[fff]]$overall 
  hold_grp <- point_estimates[[fff]]$groups
  
  if(effect_type == 'contrast'){   
    if(marginal == TRUE){
      xx <- (hold_grp[ , a1] - hold_oal[a1]) - (hold_grp[, a2] - hold_oal[a2])
    } else {
      xx <- (hold_grp[ , a1, t1] - hold_oal[a1, t1]) - (hold_grp[, a2, t2] - hold_oal[a2, t2])
    }
  } 
  else if(effect_type == 'outcome'){
    if(marginal == TRUE){
      xx <- hold_grp[ , a1] - hold_oal[a1]
    } else {
      xx <- hold_grp[  , a1, t1] - hold_oal[a1, t1]
    }
  }
  
  ## Create the V matrix for each group ##
  V_i <- array(dim = c(p+1, p+1, N))
  
  for(ii in 1:N){
    hold <- c(scores[ii, ], xx[ii])
    V_i[ , , ii] <- hold %*% t(hold)
  }
  
  ## Create final V matrix by taking the mean of each element across groups ##
  V <- apply(V_i, 1:2, sum, na.rm = T)/N

  return(V)
}
