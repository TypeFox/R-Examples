ss.aipe.cv <- function(C.of.V=NULL, width=NULL, conf.level=.95, degree.of.certainty=NULL, assurance=NULL, certainty=NULL, mu=NULL, sigma=NULL, alpha.lower=NULL, alpha.upper=NULL, Suppress.Statement=TRUE, sup.int.warns=TRUE, ...)
{
if(!is.null(certainty)& is.null(degree.of.certainty)&is.null(assurance)) degree.of.certainty<-certainty
if (is.null(assurance) && !is.null (degree.of.certainty)& is.null(certainty)) assurance <-degree.of.certainty
if (!is.null(assurance) && is.null (degree.of.certainty)& is.null(certainty)) assurance -> degree.of.certainty

if(!is.null(assurance) && !is.null (degree.of.certainty) && assurance!=degree.of.certainty) 
stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")

if(!is.null(assurance) && !is.null (certainty) && assurance!=certainty) 
stop("The arguments 'assurance' and 'certainty' must have the same value.")

if(!is.null(degree.of.certainty) && !is.null (certainty) && degree.of.certainty!=certainty) 
stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")

if(sup.int.warns==TRUE) options(warn=-1)

if(is.null(conf.level))
{
if(alpha.lower>=1 | alpha.lower<0) stop("\'alpha.lower\' is not correctly specified.")
if(alpha.upper>=1 | alpha.upper<0) stop("\'alpha.upper\' is not correctly specified.")
}

if(is.null(width)) stop("A value for \'width\' must be specified.")

if(!is.null(conf.level))
{
if(!is.null(alpha.lower) | !is.null(alpha.upper)) stop("Since \'conf.level\' is specified, \'alpha.lower\' and \'alpha.upper\' should be \'NULL\'.")
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

if(!is.null(degree.of.certainty))
{
if((degree.of.certainty <= 0) | (degree.of.certainty >= 1)) stop("The 'degree.of.certainty' must either be NULL or some value greater than .50 and less than 1.", call.=FALSE)
if(degree.of.certainty <= .50) stop("The 'degree.of.certainty' should be > .5 (but less than 1).", call.=FALSE)
}

#####

minimal.N <- 4

Lim.0 <- ci.cv(cv=cv(C.of.V=C.of.V, N=minimal.N, unbiased=TRUE), n=minimal.N, alpha.lower=alpha.lower, alpha.upper=alpha.upper, conf.level=NULL)
Current.Width <- Lim.0$Upper - Lim.0$Lower
dif <- Current.Width - width
N.0 <- minimal.N

while(dif > 0)
{
N <- N.0+1
CI.for.CV <- ci.cv(cv=cv(C.of.V=C.of.V, N=N, unbiased=TRUE), n=N, alpha.lower=alpha.lower, alpha.upper=alpha.upper, conf.level=NULL)
Current.Width <- CI.for.CV$Upper - CI.for.CV$Lower
dif <- Current.Width - width
N.0 <- N
}

#############################################################################################

# Now, incorporate a desired degree of certainty.
if(!is.null(degree.of.certainty))
{
# Here the noncentral value that will be exceeded only (1-degree.of.certainty)100% of the time is found (which leads to confidence intervals wider than desired).
beyond.CV.NCP <- qt(p=1-degree.of.certainty, df=N-1, ncp = sqrt(N)/C.of.V, lower.tail = TRUE, log.p = FALSE)

# Now the noncentral parameter is transformed into a coefficient of variation.
Lim.for.Certainty <- sqrt(N)/beyond.CV.NCP

# Now calculate sample size using the value not to be exceeded more than (1-degree.of.certainty)100% of the time.
N.gamma <- ss.aipe.cv(C.of.V=cv(C.of.V=Lim.for.Certainty, N=N, unbiased=TRUE), width=width, alpha.lower=alpha.lower, alpha.upper=alpha.upper, conf.level=NULL, degree.of.certainty=NULL, Suppress.Statement=TRUE)
}

if(is.null(degree.of.certainty)) 
{
if(Suppress.Statement==FALSE) print(paste("In order the the expected confidence interval width to be no larger than", width, ",the sample size that should be used is:", N))
return(N)
}

if(!is.null(degree.of.certainty)) 
{
if(Suppress.Statement==FALSE) print(paste("In order the the confidence interval width to be no less than", width, "with no less than", degree.of.certainty*100, "certainty, the sample size that should be used is:", N.gamma))
return(N.gamma)
}
if(sup.int.warns==TRUE) options(warn=1)
}
