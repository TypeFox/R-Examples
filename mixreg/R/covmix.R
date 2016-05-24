covmix <- function(object,x,y) {
#
# Function covmix.  To calculate the covariance matrix of the
# parameter estimates produced by mixreg().
#

theta <- object$theta
K <- length(theta)
intercept <- object$intercept
eq.var <- if(is.null(object$eq.var)) FALSE else object$eq.var
x <- as.matrix(x)
bnms <- dimnames(x)[[2]]
if(is.null(bnms)) bnms <- paste('beta',1:ncol(x),sep='')
if(intercept) {
	x <- cbind(1,x)
	bnms <- c('Int',bnms)
}
dimsok <- all(unlist(lapply(theta,function(x){length(x$beta)}))==ncol(x))
if(!dimsok) {
	cat('The values for beta are of wrong length for\n')
	cat('the dimension of the predictors.\n')
	stop('Bailing out.')
}

g      <- gfun(x,y,theta)$gamma
info.1 <- info1(x,y,theta,g)
info.2 <- info2(x,y,theta,g)

nms   <- c(outer(c(bnms,'sigsq','lambda'),1:K,paste,sep='.'))
nms   <- nms[-length(nms)]
finfo <- info.1-info.2
if(eq.var) {
	p     <- length(bnms) + 2
	ind   <- (1:K)*p - 1
	nms   <- c(nms[-ind],'sigsq')
	finfo <- aux3(finfo,ind)
}

covmat <- solve(finfo)
dimnames(covmat) <- list(nms,nms)
covmat
}
