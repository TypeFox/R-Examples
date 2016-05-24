"Rsquare2Lambda" <-
function(R2=NULL, N=NULL)
{
if(is.null(R2) | is.null(N)) stop("You must specify \'R2\' (i.e., R Square) and \'N\' (i.e., sample size) in order to calculate the noncentraility paramter.")
return((R2/(1-R2))*N)
}
