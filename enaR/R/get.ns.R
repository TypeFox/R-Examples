#' get.ns.R
#' Input = network model
#' Output = a vector of global network statistics from ena
#'
#' Borrett | July 4, 2012
#' -----------------------------------

get.ns <- function(x,balance.override=FALSE){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }

  # runs selected ena analyses that return global network statistics
  st <- enaStructure(x)$ns
  Flow <- enaFlow(x)$ns
  asc <- enaAscendency(x)
  s <- enaStorage(x)$ns
  u.f <- enaUtility(x,type='flow',eigen.check=FALSE)$ns
  u.s <- enaUtility(x,type='storage',eigen.check=FALSE)$ns

  ns <- data.frame(st,Flow,asc,s,u.f,u.s)

  return(ns)
}
  
