#-----------------------------------------------------------------------------#
# Get arguments from a function
# 
# Extracts the names of the arguments from a function, and creates a list 
# of those arguments where they exist in ... . 
# 
# @param FUN function for which to find arguments
# @param args_list a list of arguments. Defaults to NULL.
# @param ... any arguments. Those necessary for FUN must be named as appropriate for FUN
# @return list of arguments for FUN
# @export
# @examples
# myargs <- get_args(lm, formula = Sepal.Length ~ Sepal.Width, data = iris )
# summary(do.call('lm', myargs))
#-----------------------------------------------------------------------------#

get_args <- function(FUN, args_list = NULL, ...){
  dots <- append(args_list, list(...))
  arg_names <- names(formals(match.fun(FUN)))
  
  args <- dots[arg_names]
  args[sapply(args, is.null)] <- NULL
  
  return(args)
}

#-----------------------------------------------------------------------------#
# Calculate outcome mean per group per treatment level
# 
# @param Y vector of outcomes
# @param G vector of group assignments
# @param A vector of treatment assignments
# @param a value of treatment level, defaults to NA. NA is used for marginal
# estimates.
# @return matrix of group means
#-----------------------------------------------------------------------------#

group_means <- function(Y, A, G, a = NA){
  
  N <- length(unique(G))
  YA <- cbind(Y, A)
  
  vals <- by(YA, INDICES = G, function(x){
    n <- length(x[ , 1])
    
    if(is.na(a)){
      sum(x[ , 1])/n
    } else {
      sum(x[ , 1] * (x[ , 2] == a) * 1)/n
    }
  })
  
  out <- matrix(unlist(vals), nrow = N)
  
  return(out)
}

#-----------------------------------------------------------------------------#
# Create a dataset of arguments for effect estimates
#
# @param allocations vector of allocations 
# @param treatments vector of treatments. defaults to \code{c(0 ,1)}
# @return data.frame with arguments necessary for \code{\link{ipw_effect_calc}} to 
# compute all outcome, direct, indirect, total, and overall effect estimates from
# an object created from \code{\link{ipw_interference}} 
# @export
# @examples 
# effect_grid(seq(0,1, by = .1), c(0,1))
# 
#-----------------------------------------------------------------------------#

effect_grid <- function(allocations, treatments = c(0,1))
{
  marginal    <- c('TRUE', 'FALSE')
  
  # Outcomes
  g1.1 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = NA, trt2 = NA,
                    marginal = FALSE, 
                    effect_type = 'outcome', effect = 'outcome',
                    stringsAsFactors = FALSE)
  g1.2 <- expand.grid(alpha1 = allocations, trt1 = NA, 
                      alpha2 = NA, trt2 = NA,
                      marginal = TRUE, 
                      effect_type = 'outcome', effect = 'outcome',
                      stringsAsFactors = FALSE)
  
  # Direct Effects
  g2 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = NA, trt2 = treatments,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'direct',
                    stringsAsFactors = FALSE)
  g2$alpha2 <- g2$alpha1
  g2 <- g2[g2$trt1 != g2$trt2, ]
  
  # Indirect Effects
  g3 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = allocations, trt2 = NA,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'indirect',
                    stringsAsFactors = FALSE)
  g3$trt2 <- g3$trt1
  
  # Total Effects
  g4 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = allocations, trt2 = treatments,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'total',
                    stringsAsFactors = FALSE)
  g4 <- g4[g4$trt1 != g4$trt2, ]
  
  # Overall Effects
  g5 <- expand.grid(alpha1 = allocations, trt1 = NA, 
                    alpha2 = allocations, trt2 = NA,
                    marginal = TRUE, 
                    effect_type = 'contrast', effect = 'overall',
                    stringsAsFactors = FALSE)
  
  out <- rbind(g1.1, g1.2, g2, g3, g4, g5)
  rownames(out) <- NULL # Rownames aren't useful
  return(out) 
}

#-----------------------------------------------------------------------------#
#' Retreive Direct Effect estimates
#'  
#' @description Retrieves the population average direct causal effect for a 
#' specified allocation:
#' \eqn{\hat{Y}(0, alpha) - \hat{Y}(1, alpha)}{Yhat(0, alpha) - Yhat(1, alpha)}.
#'  
#' @param object an object of class \code{interference}
#' @param allocation the allocation scheme for which to estimate direct effects.
#' If NULL, then returns all direct effects.
#' @param trt.lvl1 Defaults to 0.
#' @return a data.frame with requested values
#' @export
#-----------------------------------------------------------------------------#

direct_effect <- function(object, 
                          allocation = NULL, 
                          trt.lvl1 = 0)
{
  ce  <- object$estimates[object$estimates$effect == 'direct', ]
  if(is.null(allocation)){
    out <- ce[ce$trt1 == trt.lvl1 & ce$trt2 == !trt.lvl1*1, ]
  } else{
    out <- ce[ce$alpha1 == allocation & ce$trt1 == trt.lvl1 &
              ce$alpha2 == allocation & ce$trt2 == !trt.lvl1*1, ]
  }
  rownames(out) <- NULL
  return(out)
}

#-----------------------------------------------------------------------------#
#' Retreive Indirect Effect estimates
#' 
#' @name indirect_effect
#' @description Retrieves the population average indirect causal effect for
#' specified allocations:
#' \eqn{\hat{Y}(0, alpha1) - \hat{Y}(0, alpha2)}{Yhat(0, alpha1) - Yhat(0, alpha2)}. 
#' This is the effect due to the coverage (allocation) levels.
#'  
#' @param allocation1 the allocation scheme for which to estimate indirect effects
#' @param allocation2 the allocation scheme for which to estimate indirect effects.
#' If NULL, then returns all indirect effects compared to allocation1.
#' @param trt.lvl Defaults to 0.
#' @inheritParams direct_effect
#' @return a data.frame with requested values
#' @export
#-----------------------------------------------------------------------------#

indirect_effect <- function(object, 
                                  allocation1, 
                                  allocation2 = NULL,
                                  trt.lvl = 0)
{
  ce  <- object$estimates[object$estimates$effect ==  'indirect' & 
                          object$estimates$trt1 == trt.lvl & 
                          object$estimates$trt2 == trt.lvl, ]

  if(is.null(allocation2)){
    out <- ce[ce$alpha1 == allocation1, ]
  } else {
    out <- ce[ce$alpha1 == allocation1 & ce$alpha2 == allocation2, ]
  }
  
  rownames(out) <- NULL
  return(out)
}

#-----------------------------------------------------------------------------#
#' Retreive Indirect Effect estimates
#'  
#' @rdname indirect_effect
#' @export
#-----------------------------------------------------------------------------#

ie <- indirect_effect

#-----------------------------------------------------------------------------#
#' Retrieve Total Effect estimates
#'
#' @name total_effect
#' @description Retrieves the population average total causal effect for
#' specified allocations:
#' \eqn{\hat{Y}(0, alpha1) - \hat{Y}(1, alpha2)}{Yhat(0, alpha1) - Yhat(1, alpha2)}
#'  
#' @param allocation1 the allocation scheme for which to estimate total effects
#' @param allocation2 the allocation scheme for which to estimate total effects
#' If NULL, then returns all indirect effects compared to allocation1.
#' @param trt.lvl1 Defaults to 0.
#' @inheritParams direct_effect
#' @return a data.frame with requested values
#' @export
#-----------------------------------------------------------------------------#

total_effect <- function(object, 
                               allocation1, 
                               allocation2 = NULL,
                               trt.lvl1 = 0)
{
  ce  <- object$estimates[object$estimates$effect == 'total' & 
                          object$estimates$trt1 == trt.lvl1 & 
                          object$estimates$trt2 == !trt.lvl1*1, ]
  
  if(is.null(allocation2)){
    out <- ce[ce$alpha1 == allocation1, ]
  } else {
    out <- ce[ce$alpha1 == allocation1 & ce$alpha2 == allocation2, ] 
  }

  rownames(out) <- NULL
  return(out)
}

#-----------------------------------------------------------------------------#
#' @rdname total_effect
#' @export
#-----------------------------------------------------------------------------#

te <- total_effect

#-----------------------------------------------------------------------------#
#' Retrieve Overall Effect Estimates
#' 
#' @name overall_effect
#' @description Retrieves the population average overall causal effect:
#' \eqn{\hat{Y}(alpha1) - \hat{Y}(lpha2)}{Yhat(alpha1) - Yhat(alpha2)}
#' 
#' @param allocation1 the allocation scheme for which to estimate overall effects
#' @param allocation2 the allocation scheme for which to estimate overall effects
#' @inheritParams direct_effect
#' @return a data.frame with a single row with requested values
#' @export
#-----------------------------------------------------------------------------#

overall_effect <- function(object, 
                           allocation1, 
                           allocation2 = NULL)
{
  ce  <- object$estimates[object$estimates$effect == 'overall', ]
  
  if(is.null(allocation2)){
    out <- ce[ce$alpha1 == allocation1, ]
  } else {
    out <- ce[ce$alpha1 == allocation1 & ce$alpha2 == allocation2, ]
  }

  rownames(out) <- NULL
  return(out)
}

#-----------------------------------------------------------------------------#
#' @rdname overall_effect
#' @export
#-----------------------------------------------------------------------------#

oe <- overall_effect