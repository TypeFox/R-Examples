ci.smd <- function(ncp=NULL, smd=NULL, n.1=NULL, n.2=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, ...)
{

if(is.null(ncp) & is.null(smd)) stop("You must specify either the estimated noncentral parameter 'ncp' (generally the observed t-statistic) or the standardized mean difference 'smd' (as might be obtained from the 'smd' function.).", call.=FALSE)
if(length(ncp)==1 & length(smd)==1) stop("You only need to specify either 'ncp' or 'smd', not both.", call.=FALSE)
if(is.null(n.1) | is.null(n.2)) stop("You must specify sample size per group in order to determine confidence limits.", call.=FALSE)
#if((is.null(alpha.lower) & is.null(alpha.upper)) & conf.level>=1) stop("There is a problem with your confidence level.", call.=FALSE)
if(is.null(conf.level) & sum(alpha.lower,alpha.upper)>=1) stop("There is a problem with your upper and or lower confidence limits.", call.=FALSE)

df <- n.1 + n.2 - 2

if(length(ncp)==1)
{
#if(ncp==0) stop("You need not use a noncentral method since the noncentrality parameter is zero; use the critical value from the central t-distribution.", call.=FALSE)
smd <- ncp*sqrt((n.1+n.2)/(n.1*n.2))
Limits <- conf.limits.nct(ncp, df, conf.level=conf.level, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol, ...)
Limits.L <- Limits$Lower.Limit
Limits.U <- Limits$Upper.Limit
Lower.Conf.Limit <- Limits.L*sqrt((n.1+n.2)/(n.1*n.2))
Upper.Conf.Limit <- Limits.U*sqrt((n.1+n.2)/(n.1*n.2))
Result <- list(Lower.Conf.Limit.smd=Lower.Conf.Limit, smd=smd, Upper.Conf.Limit.smd=Upper.Conf.Limit)
return(Result)
}


if(length(smd)==1)
{
#if(smd==0) stop("You need not use a noncentral method since the effect size is zero; use the critical value from the central t-distribution.", call.=FALSE)
ncp <- smd*sqrt((n.1*n.2)/(n.1+n.2))
Limits <- conf.limits.nct(ncp, df, conf.level=conf.level, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol, ...)
Limits.L <- Limits$Lower.Limit
Limits.U <- Limits$Upper.Limit
Lower.Conf.Limit <- Limits.L*sqrt((n.1+n.2)/(n.1*n.2))
Upper.Conf.Limit <- Limits.U*sqrt((n.1+n.2)/(n.1*n.2))
Result <- list(Lower.Conf.Limit.smd=Lower.Conf.Limit, smd=smd, Upper.Conf.Limit.smd=Upper.Conf.Limit)
return(Result)
}
}
