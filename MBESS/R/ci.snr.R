"ci.snr" <-
function(F.value=NULL, df.1=NULL, df.2=NULL, N=NULL, conf.level=.95,
alpha.lower=NULL, alpha.upper=NULL, ...)
{

#######################################################################33
# Preliminary information and function set-up.

if(is.null(alpha.lower) & is.null(alpha.upper))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

if ( !is.null(alpha.lower)  &&  !is.null(alpha.upper) ) 
    {
    conf.level<- 1- alpha.upper - alpha.lower
    }


if (alpha.lower>.5 || alpha.lower<0) stop (" 'alpha.lower' must be smaller than .5 and nonnegative.")
if (conf.level>1 || conf.level<0) stop (" 'conf.level' must be larger than 0 and smaller than 1. ")
if (F.value<=0) stop (" 'F.value' must be larger than 0. ")
if (N<=0 || N<=df.1+df.2) stop ("N must be larger than df.1+df.2")

if(!is.null(alpha.lower) & is.null(alpha.upper)) stop("This is a problem with the desired confidence level ('alpha.lower' specified but not 'alpha.upper').")
if(is.null(alpha.lower) & !is.null(alpha.upper)) stop("This is a problem with the desired confidence level ('alpha.upper' specified but not 'alpha.lower').")

if(is.null(df.1) | is.null(df.2) | is.null(N)) stop("You need to specify 'df.1', 'df.2', and 'N'.")
if(is.null(F.value)) stop("You must specify the observed F-value ('F.value') from the analysis of variance.")
##########################################################################

Lims <- conf.limits.ncf(F.value = F.value, conf.level = NULL, df.1 = df.1, 
        df.2 = df.2, alpha.lower = alpha.lower, alpha.upper = alpha.upper, ...)

if (!class(Lims$Lower.Limit)=="numeric") {Lims$Lower.Limit <- 0}

print(paste("The", 1-(alpha.lower + alpha.upper), "confidennce limits for the signal to noise ratio are given as:"))
return(list(Lower.Limit.Signal.to.Noise.Ratio=Lims$Lower.Limit/N, Upper.Limit.Signal.to.Noise.Ratio=Lims$Upper.Limit/N))
}