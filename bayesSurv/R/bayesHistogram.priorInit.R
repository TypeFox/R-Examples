#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       bayesHistogram.priorInit.R          ####
####                                                 ####
#### FUNCTIONS:  bayesHistogram.priorInit            ####
#########################################################

### ======================================
### bayesHistogram.priorInit
### ======================================
## Subsection of bayesHistogram
## -> just to make it more readable
##
## 26/10/2004
## 13/01/2005: check using is.null changed to check using match and is.na
## 16/01/2005: rewritten using 'give.init.Gspline', 'give.init.y', 'give.init.r'
##
bayesHistogram.priorInit <- function(prior, init, mcmc.par, design)
{
  ### Prior/init for G-spline:
  ### ========================
  gspl <- give.init.Gspline(prior, init, mcmc.par, design$dim)
  init <- attr(gspl, "init")
  prior <- attr(gspl, "prior")
  mcmc.par <- attr(gspl, "mcmc.par")  
  
  ### Prior/init for augmented observations:
  ### ======================================
  if(length(init) == 0) ininit <- "arnost"
  else                  ininit <- names(init)    
  tmp <- match("y", ininit, nomatch=NA)
  if(is.na(tmp)) init$y <- numeric(0)
  init$y <- give.init.y(init$y, design$dim, design$y.left, design$y.right, design$status)

  ### Prior/init for allocations:
  ### ===========================
  tmp <- match("r", ininit, nomatch=NA)
  if(is.na(tmp)) init$r <- numeric(0)
  init$r <- give.init.r(init$r, init$y, design$dim, prior$K, init$gamma, init$sigma, prior$c4delta, init$intercept, init$scale)

    ## Recalculate r into the vector of length n with entries from {1,...,total_length}
  if (design$dim == 1)   rsm <- init$r + prior$K[1] + 1
  else{
    if (design$dim == 2) rsm <- (init$r[,2]+prior$K[2])*(2*prior$K[1]+1) + init$r[,1]+prior$K[1]+1
    else                 stop("Unimplemented dimension appeared in bayesHistogram.priorInit()")
  }  

  ### Put everything to the resulting object:
  ### =======================================
  toreturn <- list(Gparmi = gspl$Gparmi, Gparmd = gspl$Gparmd, iter=init$iter, specification = gspl$specification,
                   y = t(init$y), r = rsm)
  attr(toreturn, "init") <- init
  attr(toreturn, "prior") <- prior
  attr(toreturn, "mcmc.par") <- mcmc.par
  
  return(toreturn)
}
