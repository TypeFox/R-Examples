rs.surv.rsadd <- 
function (formula, newdata) 
{
   
    call <- match.call()
    Terms <- terms(formula$formula)  #to rabis, ce je model mal bl smotan - as.factor ali splines ali svasta
    Terms <- delete.response(Terms)
    newdata <- model.frame(Terms,newdata)
    n <- formula$n
    if(formula$method=="max.lik"){
   	 nvar <- length(formula$coef) - length(formula$int)+1
    	 formula$coef <- formula$coef[1:nvar]
    }
    nvar <- length(formula$coef)
    nx <- nrow(newdata)
    nt <- length(formula$times)
    temp <- list(n=formula$n,time=formula$times,call=call,type="right")
    Lambda0 <- formula$Lambda0
    Lambda0 <- matrix(Lambda0,ncol=nt,nrow=nx,byrow=TRUE)
    rate <- attr(Terms, "specials")$ratetable
    R <- as.matrix(newdata[, rate,drop=FALSE])
    rat <- attributes(formula$ratetable)$dimid
    mein <- attributes(newdata[,rate])$dimnames[[2]]
    x <- match(rat,mein)
    R <- R[,x,drop=FALSE]
    newdata <- newdata[,1:nvar,drop=FALSE]
    if(any(formula$mvalue)>0)newdata <- newdata - matrix(formula$mvalue,nrow=nx,byrow=TRUE)
    R <- data.frame(R)
    names(R) <- rat
    ebx <- exp(data.matrix(newdata)%*%as.vector(formula$coef))
    ebx <- matrix(ebx,ncol=nt,nrow=length(ebx))
    Lambdae <- Lambda0*ebx
    temp$surv <- t(exp(-Lambdae))
    temp$n.event <- rep(1,nt)
    temp$n.risk <- n+1 - cumsum(temp$n.event)
    temp$time <- formula$times
    class(temp) <- c("rs.surv.rsadd", "rs.surv","survfit")
    temp
}