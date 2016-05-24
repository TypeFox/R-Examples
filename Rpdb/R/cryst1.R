#  Create an object of class 'cryst1' containing the periodical boundary conditions
#  and the name of the space group associated to a 'pdb' object.

cryst1 <- function(...)
  UseMethod("cryst1")

cryst1.default <- function(abc, abg = c(90, 90, 90), sgroup = "P1", ...)
{
  if(missing(abc)) stop("Please provide at leat 'abc'")
  to.return <- list(abc = abc, abg = abg, sgroup = sgroup)
  
  class(to.return) <- "cryst1"
  return(to.return)
}

is.cryst1 <- function(x)
{
  to.return <- any(attr(x,which="class") == "cryst1")
  return(to.return)
}