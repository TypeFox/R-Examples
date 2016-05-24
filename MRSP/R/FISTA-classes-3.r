################################################################################
#####    Functions implementing the main classes of FISTA.                 #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 10.06.2013, 18:44                                  #####
################################################################################

## control class
setClass("fista.control",
         representation = representation(
           max.iter     = "numeric",
           rel.tol      = "numeric",
           adaptive.stepsize = "ANY",
           m            = "numeric",
           BB.stepsize  = "logical",
           BB.min       = "numeric",
           BB.max       = "numeric",
           BB.every     = "numeric",
           worsen.max   = "numeric",
           trace        = "logical"),

         prototype = list(
           max.iter = 500,
           rel.tol = 1e-6,
           adaptive.stepsize = 5,
           m = 5,
           BB.stepsize = T,
           BB.min = 1e-30,
           BB.max = 1e30,
           BB.every = 5,
           worsen.max = 20,
           trace = F),

         validity = function(object){
            if(ceiling(object@max.iter) != floor(object@max.iter) |
              object@max.iter <= 0)
             return("max.iter must be a positive integer!")

           if(object@rel.tol <= 0)
             return("rel.tol has to be positive")

           return(TRUE)
         }
)

fista.control <- function(max.iter = 500, rel.tol = 1e-6, adaptive.stepsize = 5, 
                          m = 5, BB.stepsize = T, BB.min = 1e-10, BB.max = 1e30, BB.every = 5,
                          worsen.max = min(20, round(0.2 * max.iter)), trace = F){

  ## Purpose: Options for the FISTA function
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## max.iter:   maximal number of iterations of the FISTA algorithm.
  ## rel.tol:    convergence relative tolerance; the smaller the more precise.
  ## adaptive.stepsize: this controls whether an increase of the stepsize should 
  ##             be allowed or not. False turns this off, setting it to a 
  ##             numerical value means that every adaptive.stepsize iterations
  ##             during which the stepsize remained the same, it will be checked
  ##             whether the stepsize could be increased. it is recommended to 
  ##             not turn this off as the algorithm can otherwise become slow.
  ##             // currently not used!
  ## BB.stepsize: should the spectral stepsizes of Barzilai + Borwein be used?
  ## BB.min, BB.max :    lower and upper bound in safeguarding the BB.stepsize.
  ## BB.every:   e.g., BB.every = 5 means that the BB-formula is only used every
  ##             5 iterations. currently disabled!
  ## worsen.max: maximal amount of iterations in which the objective function is
  ##             worsening before fista is aborted.
  ## trace:      should several warnings and informations be plotted or not?
  ## ----------------------------------------------------------------------

  RET <- new("fista.control",
             max.iter     = max.iter,
             rel.tol      = rel.tol,
             adaptive.stepsize = adaptive.stepsize,
             m = m,
             BB.stepsize = BB.stepsize,
             BB.min = BB.min,
             BB.max = BB.max,
             BB.every = BB.every,
             worsen.max = worsen.max,
             trace = trace)
  RET
}

