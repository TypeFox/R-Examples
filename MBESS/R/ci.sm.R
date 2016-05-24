"ci.sm" <-
function(sm=NULL, Mean=NULL, SD=NULL, ncp=NULL, N=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, ...)
{
if(!is.null(Mean) | !is.null(SD))
{
if(is.null(Mean) | is.null(SD)) stop("Since either 'Mean' or 'SD' was specified, they must both be specified.")
if(!is.null(sm) | !is.null(ncp)) stop("since you specified 'Mean' and 'SD', you should not specifiy 'sm' or 'ncp'.")


if(SD < 0) stop("The estiamted standard deviation is zero, yet it should be positive.")
SM <- Mean/SD
}

if(!is.null(sm))
{
if(!is.null(Mean) | !is.null(SD) | !is.null(ncp)) stop("You need to specify only one standardized mean.", call.=FALSE)
SM <- sm
}

if(!is.null(ncp))
{
if(!is.null(Mean) | !is.null(SD) | !is.null(sm)) stop("You need to specify only one standardized mean.", call.=FALSE)
}

if(is.null(N)) stop("You must specify sample size.", call.=FALSE)


if(is.null(alpha.lower) & is.null(alpha.upper))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

if(!is.null(alpha.lower) & is.null(alpha.upper)) stop("This is a problem with the desired confidence level ('alpha.lower' specified but not 'alpha.upper').")
if(is.null(alpha.lower) & !is.null(alpha.upper)) stop("This is a problem with the desired confidence level ('alpha.upper' specified but not 'alpha.lower').")

###################################################################################################

if(!is.null(ncp))
{
SM <- ncp/sqrt(N)
}

if(!is.null(SM))
{
ncp <- SM*sqrt(N)
}

Conf.Limits <- conf.limits.nct(ncp=ncp, df=(N-1), conf.level = NULL, alpha.lower = alpha.lower, alpha.upper = alpha.upper)

LL <- Conf.Limits$Lower.Limit/sqrt(N)
UL <- Conf.Limits$Upper.Limit/sqrt(N)

Result <- list(Lower.Conf.Limit.Standardized.Mean=LL, Standardized.Mean=SM, Upper.Conf.Limit.Standardized.Mean=UL)
print(paste("The", 1-(alpha.lower+alpha.upper),"confidence limits for the standardized mean are given as:"))
return(Result)
}
