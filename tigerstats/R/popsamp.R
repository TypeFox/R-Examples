#' @title Sampling From a Population

#' @description Instructional function, and possibly a utility function for some apps.
#' 
#' @rdname popsamp
#' @usage popsamp(n,pop,...)
#' @param n number of items to sample
#' @param pop data frame, from which to sample n rows  
#' @param \ldots other arguments passed to function
#' @return The sample, as a data frame.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' data(imagpop)
#' popsamp(10,imagpop)  #Simple random sampling (no replacement)
#' popsamp(10,imagpop,replace=TRUE)  #Random sampling with replacement
popsamp <-
  function (n, pop,...) 
  {
    if (n > nrow(pop)) stop("Sample size cannot exceed population size")
    poprows <- 1:nrow(pop)
    popcols <- 1:ncol(pop)
    rsamp <- sample(poprows, size = n,...)
    "["(pop,rsamp,popcols,drop=FALSE)
    #if pop has only one variable, the default after selection of rows is to drop
    #the class from data frame to something else.  drop=FALSE prevents this.
  }

