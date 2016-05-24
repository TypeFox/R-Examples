#' enaAll --- Conduct all ecological network analyses
#' INPUT = network object
#' OUTPUT = list of analytical output
#' 
#' M. Lau | May 2013
#' ---------------------------------------------------

enaAll <- function(x = 'network object'){
  out <- list(ascendency = enaAscendency(x),
              control = enaControl(x),
              environ = enaEnviron(x),
              flow = enaFlow(x),
              mti = enaMTI(x),
              storage = enaStorage(x),
              structure = enaStructure(x),
              utility = enaUtility(x,eigen.check=FALSE))
  return(out)
}
