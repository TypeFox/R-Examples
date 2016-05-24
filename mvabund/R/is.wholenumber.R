is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
{
      return(abs(x - round(x)) < tol)
}      
