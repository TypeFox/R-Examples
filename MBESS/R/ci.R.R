ci.R <- function(R=NULL, df.1=NULL, df.2=NULL, conf.level=.95, Random.Predictors=TRUE, 
Random.Regressors, F.value=NULL, N=NULL, K=NULL, alpha.lower=NULL, alpha.upper=NULL, ...)
{

if(!is.null(R))
{
if(R<0) stop("Your multiple correlation coefficient ('R') cannot be less than zero.")
if(R>1) stop("Your multiple correlation coefficient ('R') cannot be greater than one.")
}


if(!missing(Random.Regressors) ) {Random.Predictors<- Random.Regressors}
if(missing(Random.Regressors) ) {Random.Regressors<-Random.Predictors}

if(is.null(R))
{

if(is.null(df.1)) 
{
if(is.null(K)) stop("You need to specify 'K' or 'df.1'.")
df.1 <- K
}

if(is.null(df.2))
{
if(is.null(N)) stop("You need to specify 'N' or 'df.2'.")
if(is.null(K)) stop("You need to specify 'K' or 'df.2'.")
df.2 <- N-K-1
}

R <- sqrt(F2Rsquare(F.value=F.value, df.1=df.1, df.2=df.2))

N <- NULL
K <- NULL
}


Limits <- ci.R2(R2=R^2, df.1=df.1, df.2=df.2, conf.level=conf.level, Random.Predictors=Random.Predictors, Random.Regressors=Random.Regressors, 
F.value=NULL, N=N, p=K, alpha.lower=alpha.lower, alpha.upper=alpha.upper, ...)

return(list(
Lower.Conf.Limit.R=sqrt(Limits$Lower.Conf.Limit), 
Prob.Less.Lower=Limits$Prob.Less.Lower, 
Upper.Conf.Limit.R=sqrt(Limits$Upper.Conf.Limit), 
Prob.Greater.Upper=Limits$Prob.Greater.Upper))
}
