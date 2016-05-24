#' Create population
#' 
#' Creates a new population for the ExpDE framework
#' 
#' @param popsize population size
#' @param probpars list of named problem parameters (see \code{\link{ExpDE}}).
#' 
#' @return A matrix containing the population for the ExpDE
#' 
#' @export

create_population <- function(popsize,      # population size
                              probpars)     # list of named problem parameters
{
  #Generate population of individuals within the standardized space x \in (0,1)
  
  # get problem dimension
  prob.dim <- length(probpars$xmax)
  
  return (matrix(stats::runif(n = popsize * prob.dim), 
                 nrow = popsize))
}