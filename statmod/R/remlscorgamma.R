remlscoregamma <- function(y,X,Z,mlink="log",dlink="log",trace=FALSE,tol=1e-5,maxit=40) {
#
#  Mean-dispersion fit by REML scoring for gamma responses
#  Fit ED(mu,phi) model to y with
#  g(mu)=X%*%beta and f(phi)=Z%*%gam
#
#  Gordon Smyth, Walter and Eliza Hall Institute
#  16 Dec 2002.

n <- length(y)
X <- as.matrix(X)
if(is.null(colnames(X))) colnames(X) <- paste("X",as.character(1:ncol(X)),sep="")
Z <- as.matrix(Z)
if(is.null(colnames(Z))) colnames(Z) <- paste("Z",as.character(1:ncol(Z)),sep="")
q <- dim(Z)[2]
const <- 2*sum(log(y))

# Link functions
mli <- make.link(mlink)
dli <- make.link(dlink)

# Mean family
f <- Gamma()
f$linkfun <- mli$linkfun
f$linkinv <- mli$linkinv
f$mu.eta <- mli$mu.eta
f$valideta <- mli$valideta

# initial residuals and leverages assuming constant dispersion
fitm <- glm.fit(X,y,family=f)
mu <- fitted(fitm)
d <- 2*( (y-mu)/mu - log(y/mu) )
p <- fitm$rank

# start from constant dispersion
phi <- -1/canonic.digamma(mean(d))*n/(n-p)
phi <- rep(phi,n)
fitd <- lm.fit(Z,dli$linkfun(phi))
gam <- ifelse(is.na(fitd$coef),0,fitd$coef)
if( mean(abs(fitd$residuals))/phi[1] > 1e-12 ) {
	# intercept is not in span of Z
	phi <- drop(dli$linkinv( Z %*% gam ))
	fitm <- glm.fit(X,y,weights=1/phi,mustart=mu,family=f)
	mu <- fitted(fitm)
	d <- 2*( (y-mu)/mu - log(y/mu) )
} else
	fitm <- glm.fit(X,y,weights=1/phi,mustart=mu,family=f)
dev <- const+sum(2*(lgamma(1/phi)+(1+log(phi))/phi)+d/phi)+const+2*log(prod(abs(diag(fitm$qr$qr)[1:p])))

# reml scoring
iter <- 0
if(trace) cat("Iter =",iter,", Dev =",format(dev,digits=13)," Gamma",gam,"\n")
Q2 <- array(0,c(n,p*(p+1)/2))
repeat {
	iter <- iter+1

	# gradient matrix
	eta <- dli$linkfun(phi)
	phidot <- dli$mu.eta(eta) * Z
	Z2 <- phidot / phi / sqrt(2)

	# information matrix and leverages
	Q <- qr.qy(fitm$qr, diag(1, nrow = n, ncol = p))
	j0 <- 0
	for(k in 0:(p-1)) {
		Q2[ ,(j0+1):(j0+p-k)] <- Q[ ,1:(p-k)] * Q[ ,(k+1):p]
		j0 <- j0+p-k
	}
	if(p>1) Q2[ ,(p+1):(p*(p+1)/2)] <- sqrt(2) * Q2[ ,(p+1):(p*(p+1)/2)]
	h <- drop( Q2[ ,1:p] %*% array(1,c(p,1)) )
	Q2Z <- crossprod(Q2,Z2)
	extradisp <- 2*( trigamma(1/phi) - trigamma(1/phi/h)/h )/phi^2 - (1-h)
	info <- crossprod(Z2,(extradisp+1-2*h)*Z2) + crossprod(Q2Z)

	# score vector
	deltah <- 2*(digamma(1/h/phi)+log(h)-digamma(1/phi))
	dl <- crossprod(phidot, (d - deltah)/(2*phi^2))

	# scoring step
	R <- chol(info)
	dgam <- backsolve(R,backsolve(R,dl,transpose=TRUE))
	gam <- gam + dgam

	# evaluate modified profile likelihood
	phi <- drop(dli$linkinv( Z %*% gam ))
	fitm <- glm.fit(X,y,weights=1/phi,mustart=mu,family=f)
	mu <- fitted(fitm)
	d <- 2*( (y-mu)/mu - log(y/mu) )
	dev <- const+sum(2*(lgamma(1/phi)+(1+log(phi))/phi)+d/phi)+const+2*log(prod(abs(diag(fitm$qr$qr)[1:p])))

	# iteration output
	if(trace) cat("Iter =",iter,", Dev =",format(dev,digits=13)," Gamma",gam,"\n")

	# test for convergence
	if( crossprod(dl,dgam) < tol ) break

	# test for iteration limit
	if(iter > maxit) {
		warning("Max iterations exceeded")
		break
	}
}

# Standard errors
se.gam <- sqrt(diag(chol2inv(chol(info))))
se.beta <- sqrt(diag(chol2inv(qr.R(fitm$qr))))

list(beta=fitm$coef,se.beta=se.beta,gamma=gam,se.gam=se.gam,mu=mu,phi=phi,deviance=dev,h=h)
}
