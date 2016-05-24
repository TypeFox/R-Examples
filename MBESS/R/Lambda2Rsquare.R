"Lambda2Rsquare" <-
function(Lambda=NULL, N=NULL)
{
if(is.null(Lambda) | is.null(N)) stop("You must specify \'Lambda\' (i.e., a noncentral parameter) and \'N\' (i.e., sample size) in order to calculate the noncentraility paramter.")
return(Lambda/(Lambda + N))
}
