ci.smd.c <- function(ncp=NULL, smd.c=NULL, n.C=NULL, n.E=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, ...)
{
if(is.null(ncp) & is.null(smd.c)) stop("You must specify either the estimated noncentral parameter 'ncp' (generally the observed t-statistic) or the standardized mean difference 'smd.c' (as might be obtained from the 'smd.s' function.).", call.=FALSE)
if(length(ncp)==1 & length(smd.c)==1) stop("You only need to specify either 'ncp' or 'smd.c', not both.", call.=FALSE)
if(is.null(n.C) | is.null(n.E)) stop("You must specify sample size per group in order to determine confidence limits.", call.=FALSE)
if(!is.null(conf.level) & conf.level>=1) stop("There is a problem with your confidence level.", call.=FALSE)
if(is.null(conf.level) & sum(alpha.lower,alpha.upper)>=1) stop("There is a problem with your upper and or lower confidence limits.", call.=FALSE)
if(is.null(smd.c)) smd.c <- ncp*sqrt((n.C+n.E)/(n.C*n.C))


df <- n.C-1

if(length(ncp)==1)
{
Limits <- conf.limits.nct(ncp, df, conf.level=conf.level, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol)
Limits.L <- Limits$Lower.Limit
Limits.U <- Limits$Upper.Limit
Lower.Conf.Limit <- Limits.L*sqrt((n.C+n.E)/(n.C*n.E))
Upper.Conf.Limit <- Limits.U*sqrt((n.C+n.E)/(n.C*n.E))
Result <- list(Lower.Conf.Limit.smd.c=Lower.Conf.Limit, smd.c=smd.c, Upper.Conf.Limit.smd.c=Upper.Conf.Limit)
return(Result)
}


if(length(smd.c)==1)
{
ncp <- smd.c*sqrt((n.C*n.E)/(n.C+n.E))
Limits <- conf.limits.nct(ncp, df, conf.level=conf.level, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol)
Limits.L <- Limits$Lower.Limit
Limits.U <- Limits$Upper.Limit
Lower.Conf.Limit <- Limits.L*sqrt((n.C+n.E)/(n.C*n.E))
Upper.Conf.Limit <- Limits.U*sqrt((n.C+n.E)/(n.C*n.E))
Result <- list(Lower.Conf.Limit.smd.c=Lower.Conf.Limit, smd.c=smd.c, Upper.Conf.Limit.smd.c=Upper.Conf.Limit)
return(Result)
}
}
