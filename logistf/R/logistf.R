logistf <-
function(formula = attr(data, "formula"), data = sys.parent(), pl = TRUE, alpha = 0.05,
    control, plcontrol, firth = TRUE, init, weights, plconf=NULL, dataout=TRUE, ...)
{
    #n <- nrow(data)
#    if (is.null(weights)) weights<-rep(1,nrow(data))
   call <- match.call()
   if(missing(control)) control<-logistf.control()
   if(pl==TRUE & missing(plcontrol)) plcontrol<-logistpl.control()

    mf <- match.call(expand.dots =FALSE)
    m <- match(c("formula", "data","weights", "na.action", 
        "offset"), names(mf), 0L)
 #   mf<-model.frame(formula, data=data, weights=weights)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    y <- model.response(mf)
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix 
    k <- ncol(x)    ## Anzahl Effekte
    cov.name <- labels(x)[[2]]
    weight <- as.vector(model.weights(mf)  )
    offset <- as.vector(model.offset(mf)   )
    if (is.null(offset)) offset<-rep(0,n)
    if (is.null(weight)) weight<-rep(1,n)

    if (missing(init)) init<-rep(0,k)
    if (is.null(plconf) & pl==TRUE) plconf<-1:k

    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }

    fit.full<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=1:k, init, control=control)
    fit.null<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=int, init, control=control)
    fit <- list(coefficients = fit.full$beta, alpha = alpha, terms=colnames(x), var = fit.full$var, df = (k-int), loglik =c(fit.null$loglik, fit.full$loglik),
        iter = fit.full$iter, n = sum(weight), y = y, formula = formula(formula), call=match.call(), conv=fit.full$conv)
    names(fit$conv)<-c("LL change","max abs score","beta change")
    beta<-fit.full$beta
    covs<-fit.full$var
    pi<-fit.full$pi
    fit$firth<-firth
    fit$linear.predictors <- as.vector(x %*% beta + offset)
    fit$predict <- fit.full$pi
    fit$hat.diag <- fit.full$hat.diag
    if(firth)
        fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    vars <- diag(covs)
    fit$prob <- 1 - pchisq((beta^2/vars), 1)
    fit$method.ci <- rep("Wald",k)
    fit$ci.lower <- as.vector(beta + qnorm(alpha/2) * vars^0.5)
    fit$ci.upper <- as.vector(beta + qnorm(1 - alpha/2) * vars^0.5)
    fit$alpha<-alpha
    fit$conflev<-1-alpha
    if(pl) {
        betahist.lo<-vector(length(plconf),mode="list")
        betahist.up<-vector(length(plconf),mode="list")
        pl.conv<-matrix(0,length(plconf),4)
        dimnames(pl.conv)[[1]]<-as.list(plconf)
        dimnames(pl.conv)[[2]]<-as.list(c("lower, loglik","lower, beta", "upper, loglik", "upper, beta"))
        LL.0 <- fit.full$loglik - qchisq(1 - alpha, 1)/2
        pl.iter<-matrix(0,k,2)
#        fit$ci.lower <- fit$ci.upper <- rep(0, k)
        icount<-0
        for(i in plconf) {
            icount<-icount+1
            inter<-logistpl(x, y, beta, i, LL.0, firth, -1, offset, weight, plcontrol)
            fit$ci.lower[i] <- inter$beta
            pl.iter[i,1]<-inter$iter
            betahist.lo[[icount]]<-inter$betahist
            pl.conv.lower<-t(inter$conv)
            inter<-logistpl(x, y, beta, i, LL.0, firth, 1, offset, weight, plcontrol)
            fit$ci.upper[i] <- inter$beta
            pl.iter[i,2]<-inter$iter
            betahist.up[[icount]]<-inter$betahist
            pl.conv.upper<-t(inter$conv)
            pl.conv[icount,]<-cbind(pl.conv.lower,pl.conv.upper)
            fit.i<-logistf.fit(x,y, weight=weight, offset=offset, firth, col.fit=(1:k)[-i], control=control)
            fit$prob[i] <- 1-pchisq(2*(fit.full$loglik-fit.i$loglik),1)
            fit$method.ci[i] <- "Profile Likelihood"
        }
        fit$pl.iter<-pl.iter
        fit$betahist<-list(lower=betahist.lo, upper=betahist.up)
        fit$pl.conv<-pl.conv
    }
    names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- names(fit$
        coefficients) <- dimnames(x)[[2]]
    if(dataout) {
      fit$data<-data
      fit$weights<-weight
      }
    attr(fit, "class") <- c("logistf")
    fit
}


is.logistf<-function(object){
  if(attr(object,"class")=="logistf") return(TRUE)
  else return(FALSE)
  }
  
  
coef.logistf<-function(object,...){
  return(object$coefficients)
}

confint.logistf<-function(object,parm, level=0.95, exp=FALSE, ...){
  # in fact, level is already determined in the object and will be overwritten
  level<-object$conflev
  levstr<-paste(level*100,"%",sep="")
  cimat<-cbind(object$ci.lower,object$ci.upper)
  rownames(cimat)<-names(object$ci.lower)
  colnames(cimat)<-c(paste("Lower ",levstr,sep=""),paste("Upper ",levstr,sep=""))
  if(exp) cimat<-exp(cimat)
  return(cimat)
}

vcov.logistf<-function(object,...){
  return(object$var)
  }
  
  