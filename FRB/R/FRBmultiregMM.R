
FRBmultiregMM <- function(X,...) UseMethod("FRBmultiregMM")

FRBmultiregMM.formula <- function(formula, data=NULL, ...)
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
    res <- FRBmultiregMM.default(X, Y, int = FALSE, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1L]] <- as.name("FRBmultiregMM")
    res$call <- cl
    res$xlevels <- .getXlevels(mt, mf)
    if (!is.null(miss)) res$na.action <- miss
    return(res)
}                                                 


#FRBmultiregMM.formula <- function(formula, data, ...)
#{
#    mf <- model.frame(formula=formula, data=data)
#    X <- model.matrix(attr(mf, "terms"), data=mf)
#    Y <- model.response(mf)
#
#    z <- FRBmultiregMM.default(X, Y, int = FALSE, ...)
#    return(z)
#}                                                 

FRBmultiregMM.default <- function(X, Y, int = TRUE, R=999, conf=0.95, control=MMcontrol(...),na.action=na.omit, ...)
{
# performs multivariate regression based on multivariate MM estimates, with
# fast and robust bootstrap
#
# calls: MMest_multireg(), MMboot_multireg(), MMeinfs_multireg()
#
# INPUT :
# 	Y : response matrix (n x q)
# 	X : covariates matrix (n x p) or (n x (p-1))
#   int : logical; if TRUE, an intercept column is added
#   R : number of bootstrap samples
#   conf : confidence level for bootstrap intervals
# OUTPUT :
#   res$est : (list) result of MMest_multireg()
#   res$bootest : (list) result of MMboot_multireg()
#   res$Beta : (p x q) MM estimate of the regression coefficient matrix
#   res$Sigma : (q x q) MM estimate of the error covariance matrix
#   res$SE : (p x q) bootstrap standard errors for MM-estimate Beta
#   res$CI.bca.lower : (p x q) 95% BCa lower limits for Beta
#   res$CI.bca.upper : (p x q) 95% BCa upper limits for Beta
#   res$CI.basic.lower : (p x q) 95% basic bootstrap lower limits for Beta
#   res$CI.basic.upper : (p x q) 95% basic bootstrap upper limits for Beta

# --------------------------------------------------------------------

vecop <- function(mat) {
# performs vec-operation (stacks colums of a matrix into column-vector)

nr <- nrow(mat)
nc <- ncol(mat)

vecmat <- rep(0,nr*nc)
for (col in 1:nc) {
    startindex <- (col-1)*nr+1
    vecmat[startindex:(startindex+nr-1)] <- mat[,col]
}
return(vecmat)
}

# --------------------------------------------------------------------

reconvec <- function(vec,ncol) {
# reconstructs vecop'd matrix

lcol <- length(vec)/ncol
rec <- matrix(0,lcol,ncol)
for (i in 1:ncol)
    rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]

return(rec)
}

# --------------------------------------------------------------------
# -                        main function                             -
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
q <- ncol(Y)
p <- ncol(X)
#bdp <- .5
#print(Y)

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
dimens <- p*q + q*q
MMests <- MMest_multireg(X, Y, int=FALSE, control=control)
MMBeta <- MMests$coefficients
MMSigma <- MMests$Sigma
#MMresid <- MMests$residuals  
#MMfit <- MMests$fitted.values

if (R<2) warning("argument R should be at least 2 to perform bootstrap inference; FRB is now skipped")

if (R>1) {
  bootres <- MMboot_multireg(X, Y, R=R, conf=conf, ests=MMests)

  stdsBeta <- reconvec(bootres$SE[1:(p*q)],q)
  covBeta <- bootres$cov[1:(p*q),1:(p*q)]
  lowerlimitsBeta.bca <- reconvec(bootres$CI.bca[1:(p*q),1], q)
  upperlimitsBeta.bca <- reconvec(bootres$CI.bca[1:(p*q),2], q)
  lowerlimitsBeta.basic <- reconvec(bootres$CI.basic[1:(p*q),1], q)
  upperlimitsBeta.basic <- reconvec(bootres$CI.basic[1:(p*q),2], q)
  pBeta.bca <- reconvec(bootres$p.bca[1:(p*q)], q)
  pBeta.basic <- reconvec(bootres$p.basic[1:(p*q)], q)
  if (bootres$ROK<2) {
    bootres <- NULL
    stdsBeta <- NULL
    covBeta <- NULL 
    lowerlimitsBeta.bca <- NULL
    upperlimitsBeta.bca <- NULL
    lowerlimitsBeta.basic <- NULL
    upperlimitsBeta.basic <- NULL
    pBeta.bca <- NULL
    pBeta.basic <- NULL
  }
}
else {
  bootres <- NULL
  stdsBeta <- NULL
  covBeta <- NULL 
  lowerlimitsBeta.bca <- NULL
  upperlimitsBeta.bca <- NULL
  lowerlimitsBeta.basic <- NULL
  upperlimitsBeta.basic <- NULL
  pBeta.bca <- NULL
  pBeta.basic <- NULL
}

###################################################################################

#method <- paste("Multivariate regression based on multivariate MM-estimates (bdp = ", control$bdp, ", eff = ", control$eff, ")", sep="")


#if(ncol(MMests$residuals)==1) MMresiduals=t(MMests$residuals)
#else MMresiduals=MMests$residuals

z <- list(
#Beta=MMBeta
coefficients=MMBeta, 
residuals=MMests$residuals,
fitted.values=MMests$fitted.values,
initial.coefficients=MMests$SBeta,
scale=MMests$scale,
Sigma=MMSigma, SE=stdsBeta, cov=covBeta, 
weights=MMests$weights,est=MMests, bootest=bootres, 
CI.bca.lower=lowerlimitsBeta.bca, 
        CI.bca.upper=upperlimitsBeta.bca, CI.basic.lower=lowerlimitsBeta.basic, CI.basic.upper=upperlimitsBeta.basic,
        p.bca=pBeta.bca, p.basic=pBeta.basic, conf=conf, method=MMests$method,
control=control, X=X, Y=Y, ROK=bootres$ROK,outFlag=MMests$outFlag,
df=MMests$df)
  
class(z) <- c("FRBmultireg",if (ncol(Y)==1) "lmrob") 
return(z)      

}