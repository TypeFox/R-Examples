
FRBmultiregGS <- function(X,...) UseMethod("FRBmultiregGS")

FRBmultiregGS.formula <- function(formula, data=NULL, ...)
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
    res <- FRBmultiregGS.default(X, Y, int=FALSE,...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1L]] <- as.name("FRBmultiregGS")
    res$call <- cl
    if (!is.null(miss)) res$na.action <- miss
    return(res)
}                                                 

#FRBmultiregGS.formula <- function(formula, data, ...)
#{
#    mf <- model.frame(formula=formula, data=data)
#    X <- model.matrix(attr(mf, "terms"), data=mf)
#    Y <- model.response(mf)
#
#    z <- FRBmultiregGS.default(X, Y, ...)
#    return(z)
#}                                                 

FRBmultiregGS.default <- function(X, Y, int = TRUE, R=999, bdp=0.5, conf=0.95, control=GScontrol(...),na.action=na.omit, ...)
{
# performs multivariate regression based on multivariate GS estimates, with
# fast and robust bootstrap
#
# calls: GSest_multireg.R(), GSboot_multireg.R(), GSeinfs_multireg()
#
# INPUT :
# 	Y : response matrix (n x q)
# 	X : covariates matrix (n x p), possibly including intercept column
# 		(X = as.matrix(rep(1,n)) in case of location estimation)
#   R : number of bootstrap samples
#   bdp : breakdown point of GS-estimate (determines tuning parameters)
#   conf : confidence level for bootstrap intervals
# OUTPUT :
#   res$GSest : (list) result of GSmulti.R
#   res$GSboot : (list) result of multiGSregboot.R
#   res$Beta : (p*q) GS-Beta estimate
#   res$Sigma : (q*q) GS-Sigma estimate
#   res$SE : (p*q) bootstrap standard errors for each element of GS-Beta
#   res$CI.bca.lower : (p*q) lower bounds of 95% BCa intervals for each 
#                        element of GS-Beta
#   res$CI.bca.upper : (p*q) upper bounds of 95% BCa intervals for each 
#                        element of GS-Beta
#   res$CI.basic.lower : (p*q) lower bounds of 95% basic bootstrap intervals
#                          for each element of GS-Beta
#   res$CI.basic.upper : (p*q) upper bounds of 95% basic bootstrap intervals
#                          for each element of GS-Beta
                                                           

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
Y<-YX[,1:q,drop=FALSE]
X<-YX[,-(1:q),drop=FALSE]
n <- nrow(Y)
q <- ncol(Y)
p <- ncol(X)
#bdp <- .5

if (p < 1L) stop("at least one predictor needed")
if (q < 1L) stop("at least one response needed")
if (n < (p+q)) stop("For robust multivariate regression the number of observations cannot be smaller than the 
total number of variables")

interceptdetection <- apply(X==1, 2, all)
if (any(interceptdetection)) int=TRUE
#interceptind <- (1:p)[interceptdetection==TRUE]
#colnames(X)[interceptind] <- "(intercept)"
zonderint <- (1:p)[interceptdetection==FALSE]
Xzonderint <- X[,zonderint,drop=FALSE]
X <- Xzonderint
p<-ncol(X)
if (p < 1L) stop("at least one predictor needed")

if (is.null(ynam))
    colnames(Y) <- paste("Y",1:q,sep="")
if (is.null(xnam)) 
    colnames(X) <- paste("X",1:p,sep="")

dimens <- p*q + q*q


GSests <- GSest_multireg(X, Y, int=int, bdp=bdp, control=control)
if(int) X <- cbind(rep(1,n),X)
#	GSBetawith <- GSests$coefficients
#	GSint <- t(as.matrix(GSBetawith[1,]))
#	GSBeta <- as.matrix(GSBetawith[2:(p+1),,drop=F])
	
#else{GSBeta <- GSests$coefficients}
GSBeta <- GSests$coefficients
GSSigma <- GSests$Sigma
p<-ncol(X)


if (q==1)
  {colnames(GSBeta) <- colnames(Y)
#   if(int) colnames(GSBetawith) <- colnames(Y)
   colnames(GSSigma) <- colnames(Y)
   rownames(GSSigma) <- colnames(Y)
}

if (R<2) warning("argument R should be at least 2 to perform bootstrap inference; FRB is now skipped")

if (R>1) {
#print(GSests$coef)
#print(colnames(X))
#  bootres <- GSboot_multireg(X, Y, R=R,conf=conf, ests=GSests)
  bootres <- GSboot_multireg(X, Y, R=R,conf=conf, ests=GSests) 
  stdsBeta <- reconvec(bootres$SE[1:(p*q)],q)
  covBeta <- bootres$cov[1:(p*q),1:(p*q)]
  #gives lower and upper bounds of BCA and basis bootstrap confidence intervals
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

####################################################################################

#method <- paste("Multivariate regression based on multivariate GS-estimates (breakdown point = ", bdp, ")", sep="")
#method <- list(est="GS", bdp=bdp)

z <- list(#Beta=GSBeta,
coefficients=GSBeta,
#intercept=GSint,
residuals=GSests$residuals,
fitted.values=GSests$fitted.values,
scale=GSests$scale,
Sigma=GSSigma, SE=stdsBeta,cov=covBeta,  
weights=GSests$weights, est=GSests,bootest=bootres, 
       CI.bca.lower=lowerlimitsBeta.bca,CI.bca.upper=upperlimitsBeta.bca,CI.basic.lower=lowerlimitsBeta.basic, CI.basic.upper=upperlimitsBeta.basic,
       p.bca=pBeta.bca, p.basic=pBeta.basic, conf=conf, method=GSests$method, 
control=control, X=GSests$X, Y=Y, ROK=bootres$ROK, outFlag=GSests$outFlag,
df=GSests$df)

class(z) <- "FRBmultireg"
   
return(z)

}



