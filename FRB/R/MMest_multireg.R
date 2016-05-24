MMest_multireg  <- function(X,...) UseMethod("MMest_multireg")

MMest_multireg.formula <- function(formula, data=NULL, ...)
{

# --------------------------------------------------------------------

# Returns response of formula in nice way

model.multiregresp<-function (data, type = "any") 
{
    if (attr(attr(data, "terms"), "response")) {
        if (is.list(data) | is.data.frame(data)) {
  		v <- data[[1L]]
		if (is.data.frame(data) && is.vector(v)) v <- data[,1L,drop=FALSE]
            if (type == "numeric" && is.factor(v)) {
                warning("using type=\"numeric\" with a factor response will be ignored")
            }
            else if (type == "numeric" | type == "double") 
                storage.mode(v) <- "double"
            else if (type != "any") 
                stop("invalid response type")
            if (is.matrix(v) && ncol(v) == 1L){ 
                if (is.data.frame(data)) {v=data[,1L,drop=FALSE]}
	          else {dim(v) <- NULL}}
            rows <- attr(data, "row.names")
            if (nrows <- length(rows)) {
                if (length(v) == nrows) 
                  names(v) <- rows
                else if (length(dd <- dim(v)) == 2L) 
                  if (dd[1L] == nrows && !length((dn <- dimnames(v))[[1L]])) 
                    dimnames(v) <- list(rows, dn[[2L]])
            }
            return(v)
        }
        else stop("invalid 'data' argument")
    }
    else return(NULL)
}


    mt <- terms(formula, data = data)
    if (attr(mt, "response") == 0L) stop("response is missing in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    miss <- attr(mf,"na.action")
    Y <- model.multiregresp(mf)
    Terms <- attr(mf, "terms")
    X <- model.matrix(Terms, mf)
    res <- MMest_multireg.default(X, Y, int = FALSE, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1L]] <- as.name("MMest_multireg")
    res$call <- cl
    res$xlevels <- .getXlevels(mt, mf)
    if (!is.null(miss)) res$na.action <- miss
    return(res)
}                                                 





MMest_multireg.default <- function(X, Y, int = TRUE, control=MMcontrol(...),
na.action=na.omit, ...)
{
# computes multivariate regression (M)M-estimates with auxiliary S-scale 
# INPUT:
# 	Y : response matrix (n x m)
# 	X : covariates matrix (n x p), possibly including intercept column
# 		(X = as.matrix(rep(1,n)) in case of location estimation)
# OUTPUT:
# 	res$Beta : MM-estimate of regression coefficients
#	  res$Gamma : MM-estimate of shape matrix 
#	  res$Sigma : MM-estimate of covariance matrix
# 	res$SBeta : S-estimate of regression coefficients
#   res$SGamma : S-estimate of shape matrix
#	  res$SSigma : S-estimate of covariance matrix
#   res$scale : S-estimate of scale
#
# breakdown point (bdp) set at .5
# efficiency : 95% efficiency for regression coefficients

# calls: Sest_multireg()
# --------------------------------------------------------------------
# --------------------------------------------------------------------

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes scaled Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------
# (taken from Claudia Becker's Sstart0 program)

"chi.int" <- function(p, a, c1)
return(exp(lgamma((p + a)       
        #   partial expectation d in (0,c1) of d^a under chi-squared p
  /2) - lgamma(p/2)) * 2^{a/2} * pchisq(c1^2, p + a))

# --------------------------------------------------------------------

"loceff.bw" <- function(p, c1)
{  
# called by csolve.bw.MM(); computes efficiency corresponding to c1

alpha1 <- 1/p * (chi.int(p,2,c1) - 4*chi.int(p,4,c1)/(c1^2) + 6*chi.int(p,6,c1)/(c1^4) - 4*chi.int(p,8,c1)/(c1^6) + chi.int(p,10,c1)/(c1^8))   
beta1.1 <- chi.int(p,0,c1) - 2*chi.int(p,2,c1)/(c1^2) + chi.int(p,4,c1)/(c1^4)
beta1.2 <- chi.int(p,0,c1) - 6*chi.int(p,2,c1)/(c1^2) + 5*chi.int(p,4,c1)/(c1^4)
beta1 <- (1-1/p)*beta1.1 + 1/p*beta1.2

return( beta1^2 / alpha1 )

}

# --------------------------------------------------------------------

csolve.bw.MM <- function(p, eff)
{
# constant for second Tukey Biweight rho-function for MM, for fixed location-efficiency 

maxit <- 1000
eps <- 10^(-8)
diff <- 10^6
#ctest <- csolve.bw.asymp(p,.5)
ctest <- -.4024 + 2.2539 * sqrt(p) # very precise approximation for c corresponding to 50% bdp
iter <- 1
while ((diff>eps) & (iter<maxit)) {
    cold <- ctest
    ctest <- cold * eff / loceff.bw(p,cold)
    diff <- abs(cold-ctest)
    iter <- iter+1
}
return(ctest)

}

# --------------------------------------------------------------------
# -                       main function                              -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
ynam=colnames(Y)
q=ncol(na.action(Y))
if (q < 1L) stop("at least one response needed")
X <- as.matrix(X)
xnam=colnames(X)
if (nrow(Y) != nrow(X))stop("x and y must have the same number of observations")
YX=na.action(cbind(Y,X))
Y=YX[,1:q,drop=FALSE]
X=YX[,-(1:q),drop=FALSE]
n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)
q <- ncol(Y)
if (p < 1L) stop("at least one predictor needed")
if (q < 1L) stop("at least one response needed")
if (n < (p+q)) stop("For robust multivariate regression the number of observations cannot be smaller than the 
total number of variables")

interceptdetection <- apply(X==1, 2, all)
interceptind <- (1:p)[interceptdetection==TRUE]
if (!any(interceptdetection) & int){
    X <- cbind(rep(1,n),X)
    p <- p + 1    
    interceptind <-1
    if (!is.null(xnam)) colnames(X)[1] <- "(intercept)"
}

if (is.null(ynam))
    colnames(Y) <- paste("Y",1:q,sep="")
if (is.null(xnam)) {
	colnames(X) <- paste("X",1:p,sep="")
	if (interceptdetection || int){
		colnames(X)[interceptind] <- "(intercept)"
		colnames(X)[-interceptind] <- paste("X",1:(p-1),sep="")
 	}  
}
eff <- control$eff
bdp <- control$bdp
shapeEff <- control$shapeEff
fastScontrols <- control$fastScontrols
maxiter <- control$maxIt.MM
mtol <- control$convTol.MM

c1 <- csolve.bw.MM(m, eff)
# first compute 50% breakdown S-estimator
Sresult <- Sest_multireg(X, Y, int=FALSE, bdp, fastScontrols)
auxscale <- Sresult$scale
newG <- Sresult$Gamma
newBeta <- Sresult$coefficients
newR <- Y - X %*% newBeta
psres <- sqrt(mahalanobis(newR, rep(0,m), newG))
newobj <- mean(rhobiweight(psres/auxscale,c1))
origobj <- newobj

# compute M-estimate with auxilliary scale through IRLS steps, starting
# from S-estimate
iteration <- 1
oldobj <- newobj + 1
while (((oldobj - newobj) > mtol) & (iteration < maxiter)) {
    oldobj <- newobj
    w <- scaledpsibiweight(psres/auxscale,c1)
    wbig <- matrix(rep(w,p),ncol=p)	
    wX <- X * wbig	
#    newBeta <- solve(t(wX) %*% X)%*%(t(wX) %*% Y)
    newBeta <- solve(crossprod(wX, X), crossprod(wX, Y))
    newG <- cov.wt(newR, wt=w, center=FALSE)$cov
    newG <- det(newG)^(-1/m)*newG
    newR <- Y - X %*% newBeta
    psres <- sqrt(mahalanobis(newR, rep(0,m), newG))
    newobj <- mean(rhobiweight(psres/auxscale,c1))
    iteration <- iteration+1
}

if (newobj <= origobj) {
    resBeta <- newBeta
    resshape <- newG
    rescovariance <- newG*auxscale^2
} else { # isn't supposed to happen
    resBeta <- Sresult$coefficients
    resshape <- Sresult$Gamma
    rescovariance <- Sresult$Gamma*auxscale^2
}

Fit <- X%*%resBeta
Res <- Y-Fit
method <- list(est="MM", bdp=control$bdp, eff=control$eff)
psres <- sqrt(mahalanobis(Res, rep(0,m), rescovariance))
w <- scaledpsibiweight(psres,c1)
outFlag <- (psres > sqrt(qchisq(.975, m)))
if(ncol(Res)==1) Res=t(Res)
if(ncol(Fit)==1) Fit=t(Fit)

z=list(#Beta=resBeta
coefficients=resBeta, 
residuals=Res,
fitted.values=Fit,
method=method,
control=control,
Gamma=resshape, Sigma=rescovariance, SBeta=Sresult$coefficients,
        SGamma=Sresult$Gamma, SSigma=Sresult$Gamma*auxscale^2, scale=auxscale,
df=n-(m*qr(X)$rank),X=X,Y=Y,
        c1 = c1, c0 = Sresult$c, b = Sresult$b, weights=w, outFlag=outFlag)

class(z) <- c("FRBmultireg") 
return(z)       
}
