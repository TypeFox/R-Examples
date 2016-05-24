s.u <- function(s=NULL, N=NULL, X=NULL)
{
if(is.null(X)) val <- s*((gamma((N-1)/2)*sqrt((N-1)/2))/(gamma(N/2)))

if(is.null(s) | is.null(N))
{
if(!is.null(s) | !is.null(N)) stop("Since \'s\' and \'N\' were specified, do not specify \'X\'.")
if(sum(is.na(X))!=0) stop("Missing data in \'X\' is not allowed.")

s <- (var(X))^.5
N <- length(X)

val <- s*((gamma((N-1)/2)*sqrt((N-1)/2))/(gamma(N/2)))
}
return(val)
}
