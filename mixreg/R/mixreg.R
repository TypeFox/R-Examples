mixreg <- function(x,y,ncomp=2,intercept=TRUE,eq.var=FALSE,theta.start=NULL,
                   itmax=1000,eps=1e-6,verb=TRUE,digits=7,
                   max.try=5,data.name=NULL) {
#
# Function mixreg.  To fit a mixture of regression models using the
# EM algorithm.
#

# Worry about names.
x <- as.matrix(x)
bnms <- dimnames(x)[[2]]
if(is.null(bnms)) bnms <- paste('beta',1:ncol(x),sep='')
if(intercept) {
	x <- cbind(1,x)
	bnms <- c('Int',bnms)
} else x <- as.matrix(x)

if(is.null(data.name)) {
	ynm <- deparse(substitute(y))
	xnm <- deparse(substitute(x))
	data.name <- paste(ynm,xnm,sep='.on.')
}

# Get starting values; if these are not supplied they are
# basically generated at random; it is HIGHLY recommended that
# they be supplied.  I.e. a reasonable starting guess is usually
# vital for a reasonable outcome.
K         <- ncomp
theta.old <- if(is.null(theta.start))
		init.rand(x,y,K,intercept) else theta.start
dimsok <- all(unlist(lapply(theta.old,function(x){length(x$beta)}))==ncol(x))
if(!dimsok) {
	cat('Starting values for beta are of wrong length for\n')
	cat('the dimension of the predictors.\n')
	stop('Bailing out.')
}

# Sort the initial parameter list according to the first regression
# coefficient, with the largest coefficient coming first.
tmp <- matrix(unlist(theta.old),byrow=TRUE,nrow=K)
i   <- if(intercept) 2 else 1
ind <- rev(order(tmp[,i]))
theta.old <- theta.old[ind]

# Iterate:
em.step <- 0
ntry    <- 1
theta   <- list()
sigzero <- .Machine$double.eps
repeat {
	restart <- FALSE
	em.step <- em.step + 1
	gma <- gfun(x,y,theta.old)$gamma
	lma <- apply(gma,2,mean)
	if(eq.var) sigsq <- 0
	sing <- FALSE
	for(k in 1:K) {
		nzw  <- sum(gma[,k] > sigzero)
		if(nzw > ncol(x)) {
			tmp <- lm(y ~ x - 1,weights=gma[,k])
			ccc  <- coef(tmp)
			names(ccc) <- NULL
			yhat <- fitted(tmp)
			vvv  <- sum(gma[,k]*(y-yhat)**2)
			if(eq.var) {
				sigsq <- sigsq + vvv
                		theta[[k]] <- list(beta=ccc,
                                                   sigsq=NA,lambda=lma[k])
				next
			}
			else vvv  <- vvv/sum(gma[,k])
			if(vvv < sigzero) sing <- TRUE
		}
		else sing <- TRUE
		if(sing) {
			cat('Hit singularity in likelihood surface.\n')
			if(ntry <= max.try) {
				restart <- TRUE
				cat('Trying a new starting configuration.\n')
				ntry <- ntry+1
				em.step <- 0
				theta.old <- init.rand(x,y,K,intercept)
				break
			}
			else {
				cat('Too many tries; bailing out.\n')
				return(list(theta=NA,log.like=NA,
                                         intercept=intercept,nsteps=em.step,
                                         converged=FALSE,data.name=data.name))
			}
		}
		if(restart) break
                theta[[k]] <- list(beta=ccc,sigsq=vvv,lambda=lma[k])
	}
	if(restart) next
	if(eq.var) {
		sigsq <- sigsq/length(y)
		for(k in 1:K) theta[[k]]$sigsq <- sigsq
	}
	chnge   <- max(abs(unlist(theta)-unlist(theta.old)))
	if(verb) {
		cat(paste('     EM step ',em.step,':\n',sep=''))
                cat('     max abs. change in coef.: ',
                format(round(chnge,digits)),'\n',sep='')
	}
	if(chnge < eps) {
		converged <- TRUE
		break
	}
	if(em.step == itmax) {
                cat('Failed to converge in ',itmax,' EM steps.\n',sep='')
                converged  <- FALSE
                break
        }
	theta.old <- theta
}

# Wrap it up and quit:
ll     <- gfun(x,y,theta)$log.like
M      <- K*ncol(x) + K-1 + (if(eq.var) 1 else K)
aic    <- -2*ll + 2*M
parmat <- matrix(unlist(theta),byrow=TRUE,nrow=K)
dimnames(parmat) <- list(NULL,c(bnms,'sigsq','lambda'))
rslt <- list(parmat=parmat,theta=theta,log.like=ll,aic=aic,
             intercept=intercept,eq.var=eq.var,bnms=bnms,
             nsteps=em.step,converged=converged,data.name=data.name)
class(rslt) <- 'mixreg'
rslt
}
