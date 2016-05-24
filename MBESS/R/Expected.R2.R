Expected.R2 <- function(Population.R2, N, p)
{
if(!requireNamespace("gsl", quietly = TRUE)) stop("The package 'gsl' is needed; please install the package and try again.")  
  
Value <- 1 - ((N-p-1)/(N-1))*(1-Population.R2)*gsl::hyperg_2F1(1, 1, .5*(N+1), Population.R2)
Value <- max(0, Value)
Value <- min(Value, 1)
return(Value)
}
