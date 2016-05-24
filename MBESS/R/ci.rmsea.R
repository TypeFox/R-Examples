ci.rmsea <- function(rmsea, df, N, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL)
{
if(!is.null(alpha.lower)| !is.null(alpha.upper))
    {if(is.null(alpha.upper) | is.null(alpha.lower)) stop("Both 'alpha.lower' and 'alpha.upper' must be specified")
    if(sum(alpha.lower,alpha.upper)>=1) stop("There is a problem with the specified confidence limits (their sum should not be greater than 1).",call.=FALSE)
    if (alpha.lower >= 1 | alpha.lower < 0) stop("'alpha.lower' is not correctly specified.",call.=FALSE)
    if (alpha.upper >= 1 | alpha.upper < 0) stop("'alpha.upper' is not correctly specified.",call.=FALSE)
    conf.level<-NULL  
    }
             
if(!is.null(conf.level)) 
             {
             if (!is.null(alpha.lower) | !is.null(alpha.upper))
                    {
                    stop("Since conf.level has been specified (which is done by default), you cannot specifiy 'alpha.lower' and 'alpha.upper'. If you want to specify 'alpha.lower' or 'alpha.upper', set 'conf.level=NULL'.",call.=FALSE)
                    }
             alpha.lower <- (1 - conf.level)/2
             alpha.upper <- (1 - conf.level)/2
             }

if(rmsea<0) stop("Your RMSEA cannot be less than zero, as it is defined as the minimum of 0 or the typical point estimate.", call.=FALSE)
if(df<0)    stop("Your degrees of freedom (i.e., 'df') cannot be less than zero.", call.=FALSE)
if(N<0)     stop("Your sample size (i.e., 'N') cannot be less than zero.", call.=FALSE)

chi.sq.statistic <- rmsea^2*df*(N-1)+df
#print("The Chi-Square Statistic is:")
#print(chi.sq.statistic)

#So that the (unimportant) does not print.
Orig.warn <- options()$warn
    options(warn = -1)
chi.sq.conf.limits <- conf.limits.nc.chisq(Chi.Square=chi.sq.statistic, conf.level=NULL, df=df, alpha.lower=alpha.lower, alpha.upper=alpha.upper)
options(warn=Orig.warn)
#print("The noncentral chi-square confidence limits are")
#print(chi.sq.conf.limits)

# So that the confidence interval has a real lower limit.
Low.Limit <- chi.sq.conf.limits$Lower.Limit
if(is.na(Low.Limit)==TRUE) Low.Limit <- 0
if(Low.Limit != 0) Low.Limit <- max(0, Low.Limit)

# Verify probabilities.

#Prob.Greater.Upper <- pchisq(q=chi.sq.statistic, df=df, ncp=chi.sq.conf.limits$Upper.Limit)

###########################################################################################

if(Low.Limit==0) Low.RMSEA <- 0
if(Low.Limit>0) Low.RMSEA <- (Low.Limit/(df*(N-1)))^.5
if (Low.RMSEA<=0) cat("Note: The lower confidence limit is negative and thus set to 0 based on RMSEA's definition.", "\n", "\n")
if(is.na(Low.RMSEA)) {
	Low.RMSEA <- 0
	}
Low.RMSEA <- max(0, Low.RMSEA)

#chi.sq.statistic.Low <- Low.RMSEA^2*df*(N-1)+df
#if(chi.sq.statistic.Low==0) Prob.Less.Lower <- 0
#if(chi.sq.statistic.Low>0) Prob.Less.Lower <- 1-pchisq(q=chi.sq.statistic.Low, df=df, ncp=Low.Limit)

#if(round((Prob.Less.Lower + Prob.Greater.Upper), 3) != round((alpha.lower + alpha.upper), 3))
#{
#warning("The computed confidence interval does not have the same coverage as the specified confidence interval.", call. = FALSE)
#if(alpha.lower!=0 & Prob.Less.Lower==0) warning("This is the case in this situation because the confidence interval was truncated at zero (the RMSEA cannot be #negative), and so as to not have the lower limit be negative (which would be statistically absurd), the lower limit was set to zero.", call. = FALSE)
#}

output <- list(Lower.Conf.Limit=Low.RMSEA,
	RMSEA=rmsea,
	Upper.Conf.Limit=(chi.sq.conf.limits$Upper.Limit/(df*(N-1)))^.5)
return(output)
	}