Variance.R2 <- function(Population.R2, N, p)
{
if(!requireNamespace("gsl", quietly = TRUE)) stop("The package 'gsl' is needed; please install the package and try again.")  
    
return((((N-p-1)*(N-p+1))/(N^2-1))*((1-Population.R2)^2)*(gsl::hyperg_2F1(2, 2, .5*(N+3), Population.R2)) - ((Expected.R2(Population.R2=Population.R2, N=N, p=p) - 1)^2))
}
