#-----------------------------------------------------------------------------#
#' Prints a summary of an interference object 
#'
#' @param x object of class 'interference'
#' @param ... ignored
#' @method print interference
#' @export
#-----------------------------------------------------------------------------#

print.interference <- function(x, ...)
{
  
  cols <- c('alpha1', 'trt1', 'alpha2', 'trt2', 'estimate', 'std.error', 'conf.low', 'conf.high')
  est  <- x$estimates
#Not defined for glmer class:  form <- as.character(deparse(x$models$propensity_model$formula))
  form <- x$summary$formula
  allo <- x$summary$allocations
  conf <- x$summary$conf.level
  varm <- x$summary$variance_estimation
  k    <- x$summary$nallocations
  N    <- x$summary$ngroups
  
  mina <- min(allo)
  maxa <- max(allo)

  
  de1  <- est[est$effect == "direct" & est$trt1 == 0 & est$trt2 == 1 & est$alpha1 == mina, cols ]
  de2  <- est[est$effect == "direct" & est$trt1 == 0 & est$trt2 == 1 & est$alpha1 == maxa, cols ]
  de   <- rbind(de1, de2)
  
  ie  <- est[est$effect == "indirect" & est$trt1 == 0 & est$trt2 == 0 &
                est$alpha1 == mina & est$alpha2 == maxa, cols ]
  
  te  <- est[est$effect == "total" & est$trt1 == 0 & est$trt2 == 1 &
               est$alpha1 == mina & est$alpha2 == maxa, cols ]
  
  oe  <- est[est$effect == "overall" & est$alpha1 == mina & est$alpha2 == maxa, cols ]
  
  if(length(allo) > 2){
    meda <- quantile(allo, probs = .5, type = 3)
    de3  <- est[est$effect == "direct" & est$trt1 == 0 & est$trt2 == 1 & est$alpha1 == meda, cols ]
    de   <- rbind(de, de3)
    
    ie2  <- est[est$effect == "indirect" & est$trt1 == 0 & est$trt2 == 0 &
                  est$alpha1 == mina & est$alpha2 == meda, cols ]
    ie3  <- est[est$effect == "indirect" & est$trt1 == 0 & est$trt2 == 0 &
                  est$alpha1 == meda & est$alpha2 == maxa, cols ]
    ie   <- rbind(ie, ie2, ie3)
    
    te2  <- est[est$effect == "total" & est$trt1 == 0 & est$trt2 == 1 &
                 est$alpha1 == mina & est$alpha2 == meda, cols ]
    te3  <- est[est$effect == "total" & est$trt1 == 0 & est$trt2 == 1 &
                  est$alpha1 == meda & est$alpha2 == maxa, cols ]
    te   <- rbind(te, te2, te3)
    
    oe2  <- est[est$effect == "overall" & est$alpha1 == mina & est$alpha2 == meda, cols ]
    oe3  <- est[est$effect == "overall" & est$alpha1 == meda & est$alpha2 == maxa, cols ]
    oe   <- rbind(oe, oe2, oe3)
  }
  
  de <- format(de, digits = 4)
  ie <- format(ie, digits = 4)
  te <- format(te, digits = 4)
  oe <- format(oe, digits = 4)
  
  ## Output ##
  cat(" --------------------------------------------------------------------------\n", 
      "                              Model Summary                    \n",
      "--------------------------------------------------------------------------\n",
      "Formula:", form, '\n',
      "Number of groups: ", N, '\n',
      k, "allocations were used from", mina, '(min) to', maxa, '(max) \n',
 #     "Propensity model: ", form, '\n',
      "--------------------------------------------------------------------------\n",
      "                         Causal Effect Summary                            \n",
      "                        Confidence level:", conf, "                      \n",
      "                        Variance method:", varm, "       \n",
      "--------------------------------------------------------------------------\n\n",
      "Direct Effects\n")
  print(de, row.names = FALSE)
  cat('\n', "Indirect Effects\n")
  print(ie, row.names = FALSE)
  cat('\n', "Total Effects \n")
  print(te, row.names = FALSE)
  cat('\n', 'Overall Effects \n')
  print(oe, row.names = FALSE)
  cat('\n',
      "--------------------------------------------------------------------------\n")
}
