cband <- function(object,cov.mat,x,y,alpha=0.05,xlen=100,plotit=FALSE,
                  type=NULL) {
#
# Function cband.  To do the calculations to provide 100(1-alpha)%
# confidence bands and prediction bands for a mixture of
# regressions model.
#

if(dim(as.matrix(x))[2] != 1)
	stop('Can do conf. bands for 1-var. regression only.')

theta  <- object$theta
K      <- length(theta)
int    <- object$intercept
eq.var <- object$eq.var

dimsok <- all(unlist(lapply(theta,function(x){length(x$beta)}))==1+int)
if(!dimsok) {
	cat('Values for beta are of wrong length for\n')
	cat('the dimension of the predictors.\n')
	stop('Bailing out.')
}

cdim <- length(unlist(theta)) - 1
if(eq.var) cdim <- cdim - K + 1
if(cdim != nrow(cov.mat)) {
	cat('The dimension of cov.mat is incompatible with\n')
	cat('the parameter values in object.\n')
	stop('Bailing out.')
}

if(is.null(type)) type <- 'both'
alpha.use <- if(type=='both') alpha/2 else alpha

xf <- if(int) cbind(1,seq(min(x),max(x),length=xlen))
	else
		as.matrix(seq(min(x),max(x),length=xlen))
tv <- qnorm(1-alpha.use)

do.up <- switch(type,both=TRUE,upper=TRUE,lower=FALSE)
if(is.null(do.up)) {
	cat('Argument type must be one of upper, lower, or both.\n')
	stop('Bailing out.')
}
do.dn <- switch(type,both=TRUE,upper=FALSE,lower=TRUE)
bnds <- list()
for(k in 1:K) {
	beta <- theta[[k]]$beta
	if(eq.var)
		ind  <- if(int) 3*(k-1) + 1:2 else 2*(k-1) + 1
	else
		ind  <- if(int) 4*(k-1) + 1:2 else 3*(k-1) + 1
	yf   <- xf%*%beta
	vf   <- apply(xf*(xf%*%cov.mat[ind,ind]),1,sum)
	ucb  <- if(do.up) yf + tv*sqrt(vf) else NULL
	lcb  <- if(do.dn) yf - tv*sqrt(vf) else NULL
	upb  <- if(do.up) yf + tv*sqrt(theta[[k]]$sigsq + vf) else NULL
	lpb  <- if(do.dn) yf - tv*sqrt(theta[[k]]$sigsq + vf) else NULL
	bnds[[k]] <- cbind(lcb,ucb,lpb,upb)
}

if(int) xf <- xf[,2] else xf <- c(xf)
rslt <- list(theta=theta,intercept=int,x=x,y=y,xf=xf,bnds=bnds,
             type=type,alpha=alpha)
class(rslt) <- 'cband'
if(plotit) {
	plot(rslt)
	return(invisible(rslt))
} else return(rslt)
}
