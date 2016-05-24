#'VAR representation
#'
#'Show the VAR representation of a VECM
#'
#'
#'@aliases VARrep VARrep.VECM VARrep.VAR
#'@param object An object of class \sQuote{VECM} created by \code{\link{VECM}},
#'or of class \sQuote{VAR} created by \code{\link{lineVar}}
#'@param \dots Currently not used
#'@return A matrix containing the parameters of the VECM under their VAR
#'representation.
#'@author Matthieu Stigler
#'@references Hamilton (1994) \emph{Time Series Analysis}, Princeton University
#'Press
#'@keywords ts VECM VAR cointegration
#'@examples
#'
#'
#'data(barry)
#'
#'# VECM model:
#'mod_vecm <- VECM(barry, lag=2, estim="ML")
#'VARrep(mod_vecm)
#'
#'# VAR model:
#'mod_var <- lineVar(barry, lag=2, I="diff")
#'VARrep(mod_var)
#'
#'

#' @export
VARrep  <- function (object, ...)  
  UseMethod("VARrep")

#' @rdname VARrep
#' @method VARrep VECM
#' @S3method VARrep VECM
VARrep.VECM <- function(object, ...) {

  lag <- object$lag
  k <- object$k
  r <- object$model.specific$r
  co <- object$coefficients
  include <- object$include
  LRinclude <- object$model.specific$LRinclude

##obtain matrices
  betas <- object$model.specific$beta
  if(LRinclude!="none") betas <- betas[1:k,, drop=FALSE]

  ## Pi matrix
  Pi <-  co[, grep("ECT", colnames(co))]%*%t(betas)
  ## A_i matrix
  Amat <- matrix(NA, nrow=k, ncol=lag*(k)+k)

  if(lag>0){
    ## A_lag+1 matrix
    Amat[,(1:k)+k*lag] <- -co[,grep(paste("-", lag, "$", sep=""), colnames(co))]
    ## A_lag+1 matrix
    if(lag>1) for(i in 1:(lag-1)) Amat[,(1:k)+k*i] <- -(co[,grep(paste("-", i, "$", sep=""), colnames(co))] -co[,grep(paste("-", i+1, "$", sep=""), colnames(co))])

    cumulMat <- matrix(0, k,k)
    for(i in 1:lag) cumulMat <- cumulMat + Amat[,(1:k)+k*i]
    Amat[, 1:k] <- Pi + (diag(k)- cumulMat )
  } else {
    Amat[, 1:k] <- Pi + diag(k)
  }
## Names
  varNames <- colnames(object$model)[1:k]
  colnames(Amat) <- paste(rep(varNames, lag+1), rep(1:(lag+1), each=k), sep=".l")

## Add deterministic terms
  if(include!="none"){
    incName <- switch(object$include, "const"="Intercept", trend="Trend", both="Intercept|Trend")
    incVar <- co[,grep(incName , colnames(co)),drop=FALSE]
    if(LRinclude!="none"){
      Pi_all <-  co[, grep("ECT", colnames(co))]%*%t(object$model.specific$beta)
      Pi_deter <- Pi_all[,"trend", drop=FALSE]
      colnames(Pi_deter) <- "Trend"
      Amat <- cbind(Pi_deter,Amat)
    }
    Amat <- cbind(incVar,Amat)
    colnames(Amat) <- gsub("Intercept", "constant", colnames(Amat))
    
  } else if(LRinclude!="none"){
    Pi_all <-  co[, grep("ECT", colnames(co))]%*%t(object$model.specific$beta)
    Pi_deter <- Pi_all[,switch(LRinclude, "const"="const", "trend"="trend", "both"=c("const", "trend")), drop=FALSE]
    colnames(Pi_deter) <- switch(LRinclude, "const"="constant", "trend"="Trend", "both"=c("constant", "Trend"))
    Amat <- cbind(Pi_deter,Amat)
  }
  rownames(Amat) <- gsub("Equation ","",rownames(co))

## Add exogen terms
  if(object$exogen){
    co_exo <- co[,tail(1:ncol(co),object$num_exogen)]
    Amat <- cbind(Amat, co_exo)
  }

## res
  Amat
}

#' @rdname VARrep
#' @method VARrep VAR
#' @S3method VARrep VAR
VARrep.VAR <- function(object, ...) {

  I <- attr(object, "varsLevel")

  if(I=="ADF"){
    stop("Sorry, VARrep not yet implemented for type=ADF. Please use corresponding level formulation with lag+1")
  } else if(I=="level"){
    res <- coef(object)
  } else if(I=="diff"){
    lag <- object$lag
    k <- object$k
    co <- coef(object)
    include <- object$include
    origNames <- colnames(object$model[,1:k])

    comat <- matrix(NA, ncol=k*(lag+1), nrow=k)

  ## first lag
    comat[,(1:k)] <- diag(k)

    for(i in 1:lag){
      comat[,(1:k)+k*(i-1)] <- comat[,(1:k)+k*(i-1)]+
        co[,grep(paste("-", i, "$", sep=""), colnames(co))] 
#       if(i>1){
	comat[,(1:k)+k*i] <- -co[,grep(paste("-", i, "$", sep=""), colnames(co))] 
#       }
    }
  ## names
    colnames(comat) <- paste(rep(origNames, lag+1), rep(1:(lag+1), each=k), sep=".l")

  ## add deterministic terms
    if(include!="none"){
      inc_name <- switch(include, "none"=NULL, "const"="Intercept", "trend"="Trend", "both"=c("Intercept","Trend"))
      comat <- cbind(co[,inc_name,drop=FALSE], comat)
    }
    res <- comat
  ## Add exogen terms
    if(object$exogen){
      co_exo <- co[,tail(1:ncol(co),object$num_exogen)]
      res <- cbind(res, co_exo)
    }

  } 


##
return(res)
}
  


############################################################
#################### vec2var.tsDyn 
############################################################



vec2var.tsDyn <- function(x){

  model <- if(inherits(x,"VECM")) "VECM" else "VAR"
  co <- coef(x)
  lag <- ifelse(model=="VECM",x$lag+1, x$lag)
  K <- x$k
  include <- x$include

## VECM case: 
  if(model=="VECM"){
    LRinclude <- x$model.specific$LRinclude
    if(LRinclude!="none"){
      if(LRinclude=="const"){
	include <- "const"
      } else if(LRinclude=="trend"){
	include <- if(include=="const") "both" else "trend"
      } else if(LRinclude=="both"){
	include <- "both"
      }
    }
  }

## Take vec2var representation for VECMs
  if(model=="VECM") {
    co <- VARrep(x)
  }
  rownames(co) <- gsub("Equation ", "", rownames(co))
  colnames(co) <- gsub(" -([0-9]+)","\\.l\\1", colnames(co))
  colnames(co) <- gsub("Intercept","constant", colnames(co))

## detcoeffs
  detcoeffs <- co[,grep("constant|Trend", colnames(co)), drop=FALSE]

## A
  A <- list()
  for(i in 1:lag) A[[i]] <- co[,grep(paste("\\.l", i, sep=""), colnames(co)), drop=FALSE]
  names(A) <- paste("A", 1:lag, sep="")

## Rank
  rank <- if(model=="VECM") x$model.specific$r else K

## vecm
  ecdet <- if(model=="VECM") x$model.specific$LRinclude else "none"
  aChar <- "fakeChar"
  vecm<- new("ca.jo", season = NULL, dumvar=NULL, ecdet=ecdet,lag=as.integer(lag),spec="transitory", lambda=aChar)

## datamat
  if(model=="VAR"){
    datamat <- as.matrix(tail(as.data.frame(x$model),-lag))
  } else {
    newx <- lineVar(x$model[,1:K], lag=lag, include=include)
    datamat <- as.matrix(tail(as.data.frame(newx$model),-lag))
  }
  colnames(datamat) <- gsub(" -([0-9]+)","\\.l\\1", colnames(datamat))
  colnames(datamat) <- gsub("Intercept","constant", colnames(datamat))

## residuals
  resids <- residuals(x)
  colnames(resids) <- paste("resids of", colnames(resids))
## Return:
  result <- list(deterministic = detcoeffs, A = A, p = lag, K = K, y = as.matrix(x$model[,1:x$k]), obs = x$t, totobs = 
		  x$T, call = match.call(), vecm = vecm, datamat = datamat, resid = resids, r = rank)

  class(result) <- "vec2var"
  return(result)   

}


############################################################
#################### Methods
############################################################

predictOld.VAR <- function(object,...){
  if(object$include%in%c("both","none")) stop("Does not work with include='none' or 'both'")
  if(attr(object, "varsLevel")!="level") stop("Does not work with VAR in diff or ADf specification")

  predict(vec2var.tsDyn(object), ...)
}

predictOld.VECM <- function(object,...){
  if(object$include=="none"&&object$model.specific$LRinclude=="none") stop("Does not work with include='none'")
 predict(vec2var.tsDyn(object), ...)
}

#' @export irf
#' @S3method irf nlVar
irf.nlVar <- function(x, impulse=NULL, response=NULL, n.ahead=10, ortho=TRUE, cumulative=FALSE, boot=TRUE, ci=0.95, runs=100, seed=NULL, ...){
  model <- attr(x, "model")
  if(model=="VECM"){
    LRinc <- x$model.specific$LRinclude
    inc <- x$include
    if(LRinc=="both"|inc=="none"&LRinc=="none") stop("Sorry, irf() is not available for this specification of deterministic terms.")
  }
 irf(vec2var.tsDyn(x), impulse=impulse, response=response, n.ahead = n.ahead, ortho=ortho, cumulative=cumulative, boot=boot, ci=ci, runs=runs, seed=seed, ...)
}

#' @export fevd
#' @S3method fevd nlVar
fevd.nlVar <- function(x, n.ahead=10, ...){
  model <- attr(x, "model")
  if(model=="VECM"){
    LRinc <- x$model.specific$LRinclude
    inc <- x$include
    if(LRinc=="both"|inc=="none"&LRinc=="none") warning("Not guaranted to work with this specification of deterministic terms.")
  }
 fevd(vec2var.tsDyn(x),n.ahead=n.ahead, ...)
}


####### Predict 


predict.VECMMiddleold <- function(object, newdata, n.ahead=5, ...){
  lag <- object$lag
  k <- object$k
  include <- object$include
  LRinclude <- object$model.specific$LRinclude

## get VAR rrepresentation
  B <- VARrep(object)

## check deterministc specification
  if(LRinclude!="none"){
    if(LRinclude=="const"){
      include <- "const"
    } else if(LRinclude=="trend"){
      include <- if(include=="const") "both" else "trend"
    } else if(LRinclude=="both"){
      include <- "both"
    }
  }
## setup starting values (data in y), innovations (0)
  original.data <- object$model[,1:k]
  starting <-  myTail(original.data,lag+1) 
  innov <- matrix(0, nrow=n.ahead, ncol=k)

  if(!missing(newdata)) {
    if(!inherits(newdata, c("data.frame", "matrix","zoo", "ts"))) stop("Arg 'newdata' should be of class data.frame, matrix, zoo or ts")
    if(nrow(newdata)!=lag+1) stop("Please provide newdata with nrow=lag+1 (note lag=p in VECM representation corresponds to p+1 in VAR rep)")
    starting <-  newdata 
  }

## use VAR sim
  res <- VAR.sim2(B=B, lag=lag+1, n=n.ahead, starting=starting, innov=innov, include=include)

## results
  colnames(res) <- colnames(original.data )
  res <- tail(res, n.ahead)
  rownames(res) <- (nrow(original.data)+1):(nrow(original.data)+n.ahead)

  return(res)
}


myHead <- function(x, n=6){

  if(inherits(x, "ts")){
    res <-  apply(x,2,head,n)
    if(n==1) res <- matrix(res, nrow=1)
  } else {
    res <-  head(x,n) 
  }

  res
}

myTail <- function(x, n=6,...){

  if(inherits(x, "ts")){
    res <-  apply(x,2,tail,n,...)
    if(n==1) res <- matrix(res, nrow=1)
  } else {
    res <-  tail(x,n,...) 
  }

  res
}

############################################################
#################### EXAMPLES, tests
############################################################

if(FALSE){


library(tsDyn)

#data(zeroyld)
vec1 <- VECM(zeroyld, lag=2, estim="ML")
predict(vec1 )
tsDyn:::predictOld.VECM(vec1, n.ahead=5)
fevd(vec1 )
irf(vec1, runs=10 )

varpToDf <- function(x) matrix(sapply(x$fcst, function(x) x[,"fcst"]), nrow=nrow(x$fcst[[1]]))


### Comparisons
library(vars)
#data(Canada)
n <- nrow(Canada)

VECM_tsD <- VECM(Canada, lag=2, estim="ML")
VAR_tsD <- lineVar(Canada, lag=2)
VAR_tsD_tovars <-tsDyn:::vec2var.tsDyn(VAR_tsD)
VECM_tsD_tovars <-tsDyn:::vec2var.tsDyn(VECM_tsD)

VAR_vars <- VAR(Canada, p=2)
VECM_vars1 <- cajorls(ca.jo(Canada, K=3, spec="transitory"))
VECM_vars <- vec2var(ca.jo(Canada, K=3, spec="transitory"))


### Compare VECM methods:

### predict: OK!!
all.equal(predict(VECM_tsD)$fcst,predict(VECM_vars)$fcst)
all.equal(predict(VECM_tsD)$endog,predict(VECM_vars)$endog)

### fevd: OK!!
all.equal(fevd(VECM_tsD),fevd(VECM_vars))

### irf: OK!!
all.equal(irf(VECM_tsD, boot=FALSE)$irf,irf(VECM_vars, boot=FALSE)$irf)
all.equal(irf(VECM_tsD, boot=FALSE)$Lower,irf(VECM_vars, boot=FALSE)$Lower)

all.equal(irf(VECM_tsD, boot=TRUE,runs=2, seed=1234)$Lower, irf(VECM_vars, boot=TRUE,runs=2, seed=1234)$Lower)
all.equal(irf(VECM_tsD, boot=TRUE,runs=2, seed=1234)$Upper, irf(VECM_vars, boot=TRUE,runs=2, seed=1234)$Upper)


### Compare VECM methods:
predict(VAR_tsD)
all.equal(predict(VAR_tsD)$fcst,predict(VAR_vars)$fcst)
all.equal(predict(VECM_tsD)$endog,predict(VECM_vars)$endog)


### compare VARrep
#data(denmark)
dat_examp <- denmark[,2:3]


toVARrep <- function(ca.jo){
  vec2<- vec2var(ca.jo)
  lags <- vec2$A[[1]]
  if(length(vec2$A)>1) for(i in 2:length(vec2$A)) lags <- cbind(lags, vec2$A[[i]])
  cbind(vec2$deterministic,lags)
}

toVARrep(ca.jo=ca.jo(dat_examp,  K=2, spec="transitory"))

VARrep(VECM(dat_examp, lag=1, include="const", estim="ML"))
toVARrep(ca.jo(dat_examp,  K=2, spec="transitory"))

all.equal(VARrep(VECM(dat_examp, lag=1, include="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=2, spec="transitory")), check.attributes=FALSE)
all.equal(VARrep(VECM(dat_examp, lag=2, include="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=3, spec="transitory")), check.attributes=FALSE)

all.equal(VARrep(VECM(dat_examp, lag=1, LRinclude="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=2, spec="transitory", ecdet="const")), check.attributes=FALSE)
all.equal(VARrep(VECM(dat_examp, lag=2, LRinclude="const", estim="ML")), toVARrep(ca.jo(dat_examp,  K=3, spec="transitory", ecdet="const")), check.attributes=FALSE)

all.equal(VARrep(VECM(dat_examp, lag=1, LRinclude="trend", estim="ML")), toVARrep(ca.jo(dat_examp,  K=2, spec="transitory", ecdet="trend")), check.attributes=FALSE)
all.equal(VARrep(VECM(dat_examp, lag=2, LRinclude="trend", estim="ML")), toVARrep(ca.jo(dat_examp,  K=3, spec="transitory", ecdet="trend")), check.attributes=FALSE)

all.equal(VARrep(VECM(Canada, lag=1, LRinclude="trend", estim="ML")), toVARrep(ca.jo(Canada,  K=2, spec="transitory", ecdet="trend")), check.attributes=FALSE)
VECM(Canada, lag=1, LRinclude="trend", estim="ML")
cajorls(ca.jo(Canada,  K=2, spec="transitory", ecdet="trend"))$rlm



#### compare slots: VECM
all.equal(VECM_tsD_tovars,VECM_vars)
all.equal(VECM_tsD_tovars$deterministic,VECM_vars$deterministic)
all.equal(VECM_tsD_tovars$A,VECM_vars$A)
all.equal(VECM_tsD_tovars$y,VECM_vars$y)
all.equal(VECM_tsD_tovars$resid,VECM_vars$resid)
attributes(VECM_vars$resid)
attributes(VECM_tsD_tovars$resid)

all.equal(VECM_tsD_tovars$datamat,VECM_vars$datamat)
all.equal(VECM_tsD_tovars$p,VECM_vars$p)
all.equal(VECM_tsD_tovars$r,VECM_vars$r)


#### compare slots: VAR


## compare coefs
coef(VECM_tsD)
t(coef(VECM_vars1$rlm))

## compare vec2var
VECM_vars$A$A3
vec2var.tsDyn(VECM_tsD)

## compare vec2var.tsDyn
VECM_tsD_tovars$A$A3
VECM_vars$A$A3

## compare residuals
head(VECM_vars$datamat)
head(VECM_tsD_tovars$datamat)

head(residuals(VECM_tsD_tovars),2)
head(residuals(VECM_vars),2)

## compare predict
VECM_tsD_tovars$r

predict(VECM_tsD_tovars, n.ahead=5)$fcst$U
predict(VECM_vars, n.ahead=5)$fcst$U

### IRF

head(irf(VECM_tsD_tovars, boot=FALSE)$irf$U,3)
head(irf(VECM_vars, boot=FALSE)$irf$U,3)

head(irf(VECM_tsD_tovars, boot=TRUE)$Upper$U,3)
head(irf(VECM_vars, boot=TRUE)$Upper$U,3)

### FEVD
head(fevd(VECM_tsD_tovars)$U,3)
head(fevd(VECM_vars)$U,3)

### predict


### compare prediction and actual
Var_1 <- lineVar(Canada, lag=1)
all.equal(predict(Var_1),predict(Var_1, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1),1),predict(Var_1, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Var_2 <- lineVar(Canada, lag=2)
all.equal(predict(Var_2),predict(Var_2, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(tail(fitted(Var_2),1),predict(Var_2, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Var_1_t <- lineVar(Canada, lag=1, include="trend")
all.equal(predict(Var_1_t),predict(Var_1_t, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1_t),1),predict(Var_1_t, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Var_1_no <- lineVar(Canada, lag=1, include="none")
all.equal(predict(Var_1_no),predict(Var_1_no, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1_no),1),predict(Var_1_no, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Var_1_bo <- lineVar(Canada, lag=1, include="both")
all.equal(predict(Var_1_bo),predict(Var_1_bo, newdata=Canada[n,,drop=FALSE]))
all.equal(tail(fitted(Var_1_bo),1),predict(Var_1_bo, n.ahead=1, newdata=Canada[(n-1),,drop=FALSE]), check.attributes=FALSE)

Var_1_dif <- lineVar(Canada, lag=1, I="diff")
all.equal(predict(Var_1_dif),predict(Var_1_dif, newdata=Canada[(n-1):n,,drop=FALSE]))
all.equal(tail(fitted(Var_1_dif, level="original"),1),predict(Var_1_dif, n.ahead=1, newdata=Canada[(n-2):(n-1),,drop=FALSE]), check.attributes=FALSE)


### VECM
Vecm_1_co <- VECM(Canada, lag=1, include="const", estim="ML")
Vecm_1_co_vars <- ca.jo(Canada, K=2,spec="transitory",  ecdet="none")
all.equal(predict2.VECM(Vecm_1_co),predict2.VECM(Vecm_1_co, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_1_co), varpToDf( predict(Vecm_1_co,n.ahead=1)), check.attributes=FALSE)
all.equal(predict2.VECM(Vecm_1_co), varpToDf( predict(vec2var(Vecm_1_co_vars ),n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_1_co, level="original"),1),predict2.VECM(Vecm_1_co, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)
all.equal(tail(fitted(vec2var(Vecm_1_co_vars)),1),predict2.VECM(Vecm_1_co, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)


Vecm_2_co <- VECM(Canada, lag=2, include="const", estim="ML")
Vecm_2_co_vars <- ca.jo(Canada, K=3,spec="transitory",  ecdet="none")
all.equal(predict2.VECM(Vecm_2_co),predict2.VECM(Vecm_2_co, newdata=Canada[c(n-2,n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_2_co), varpToDf( predict(Vecm_2_co,n.ahead=1)), check.attributes=FALSE)
all.equal(predict2.VECM(Vecm_2_co), varpToDf( predict(vec2var(Vecm_2_co_vars ),n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_2_co, level="original"),1),predict2.VECM(Vecm_2_co, n.ahead=1, newdata=Canada[c(n-3,n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_tr <- VECM(Canada, lag=1, include="trend", estim="ML")
all.equal(predict2.VECM(Vecm_1_tr),predict2.VECM(Vecm_1_tr, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_1_tr), varpToDf( predict(Vecm_1_tr,n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_1_tr, level="original"),1),predict2.VECM(Vecm_1_tr, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_non <- VECM(Canada, lag=1, include="none", estim="ML")
all.equal(predict2.VECM(Vecm_1_non),predict2.VECM(Vecm_1_non, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_non, level="original"),1),predict2.VECM(Vecm_1_non, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)

Vecm_1_LRco <- VECM(Canada, lag=1, LRinclude="const", estim="ML")
Vecm_1_LRco_vars <- ca.jo(Canada, K=2,spec="transitory",  ecdet="const")
all.equal(predict2.VECM(Vecm_1_LRco),predict2.VECM(Vecm_1_LRco, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_1_LRco ), varpToDf( predict(Vecm_1_LRco ,n.ahead=1)), check.attributes=FALSE)
all.equal(predict2.VECM(Vecm_1_LRco), varpToDf( predict(vec2var(Vecm_1_LRco_vars),n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_1_LRco, level="original"),1),predict2.VECM(Vecm_1_LRco, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)


Vecm_1_LRt <- VECM(Canada, lag=1, LRinclude="trend", estim="ML")
Vecm_1_LRt_vars <- ca.jo(Canada, K=2,spec="transitory",  ecdet="trend")
all.equal(predict2.VECM(Vecm_1_LRt),predict2.VECM(Vecm_1_LRt, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(predict2.VECM(Vecm_1_LRt ), varpToDf( predict(Vecm_1_LRt ,n.ahead=1)), check.attributes=FALSE)
all.equal(predict2.VECM(Vecm_1_LRt), varpToDf( predict(vec2var(Vecm_1_LRt_vars),n.ahead=1)), check.attributes=FALSE)
all.equal(tail(fitted(Vecm_1_LRt, level="original"),1),predict2.VECM(Vecm_1_LRt, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)


Vecm_1_LRbo <- VECM(Canada, lag=1, LRinclude="both", estim="ML")
all.equal(predict2.VECM(Vecm_1_LRbo),predict2.VECM(Vecm_1_LRbo, newdata=Canada[c(n-1,n),,drop=FALSE]))
all.equal(tail(fitted(Vecm_1_LRbo, level="original"),1),predict2.VECM(Vecm_1_LRbo, n.ahead=1, newdata=Canada[c(n-2,n-1),,drop=FALSE]), check.attributes=FALSE)


predict(lineVar(Canada, lag=1, include="none"))
predict2.VECM(VECM(Canada, lag=1))
predict2.VECM(VECM(Canada, lag=1), newdata=Canada[(nrow(Canada)-1):nrow(Canada),,drop=FALSE])




all.equal(sapply(predict(VECM_tsD, n.ahead=5)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=5), check.attributes=FALSE, tol=1e-07)
all.equal(sapply(predict(VECM_tsD, n.ahead=10)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=10), check.attributes=FALSE, tol=1e-07)


### VAR
all.equal(predict(lineVar(Canada, lag=2), n.ahead=3), sapply(predict(VAR(Canada, p=2), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict(lineVar(Canada, lag=4), n.ahead=3), sapply(predict(VAR(Canada, p=4), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)


all.equal(predict(lineVar(Canada, lag=1, include="none"), n.ahead=3), sapply(predict(VAR(Canada, p=1, type="none"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict(lineVar(Canada, lag=2, include="none"), n.ahead=3), sapply(predict(VAR(Canada, p=2, type="none"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)

all.equal(predict(lineVar(Canada, lag=1, include="trend"), n.ahead=3), sapply(predict(VAR(Canada, p=1, type="trend"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict(lineVar(Canada, lag=2, include="trend"), n.ahead=3), sapply(predict(VAR(Canada, p=2, type="trend"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)

all.equal(predict(lineVar(Canada, lag=1, include="both"), n.ahead=3), sapply(predict(VAR(Canada, p=1, type="both"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)
all.equal(predict(lineVar(Canada, lag=2, include="both"), n.ahead=3), sapply(predict(VAR(Canada, p=2, type="both"), n.ahead=3)$fcst, function(x) x[,"fcst"]), check.attributes=FALSE)

### VECM
all.equal(sapply(predict(VECM_tsD, n.ahead=5)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=5), check.attributes=FALSE)
all.equal(sapply(predict(VECM_tsD, n.ahead=10)$fcst, function(x) x[,"fcst"]), predict2(VECM_tsD, n.ahead=10), check.attributes=FALSE)


ve2_tsD <- VECM(Canada, lag=2, estim="ML")
ve2_var <- vec2var(ca.jo(Canada, K=3, spec="transitory"))
all.equal(sapply(predict(ve2_tsD, n.ahead=5)$fcst, function(x) x[,"fcst"]), predict2(ve2_tsD, n.ahead=5), check.attributes=FALSE)
all.equal(sapply(predict(ve2_var, n.ahead=10)$fcst, function(x) x[,"fcst"]), predict2(ve2_tsD, n.ahead=10), check.attributes=FALSE)


coef(lineVar(Canada, lag=2, I="diff"))
VARrep.VAR(lineVar(Canada, lag=2, I="diff"))

coef(lineVar(Canada, lag=1, I="diff"))
VARrep.VAR(lineVar(Canada, lag=1, I="diff"))



}



