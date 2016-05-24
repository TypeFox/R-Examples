# this function is originally named BYlogreg. It is written by Christophe Croux and Gentiane Haesbroeck, also collected in robustbase package after certain version.
# to avoid dependence on rrcov and numerical instability from that, I hard code initwml to FALSE

rlogit  <- function(formula, dat, const=0.5, kmax=1e3, maxhalf=10, verbose=FALSE)
{
    initwml=FALSE
     
#  Computation of the estimator of Bianco and Yohai (1996) in logistic regression
#  -------------
#  Christophe Croux, Gentiane Haesbroeck
#  (thanks to Kristel Joossens and Valentin Todorov for improving the code)
#
#  This program computes the estimator of Bianco and Yohai in
#  logistic regression. By default, an intercept term is included
#  and p parameters are estimated.
#
#  For more details we refer to
#     Croux, C., and Haesbroeck, G. (2003), ``Implementing the Bianco and Yohai estimator for Logistic Regression'',
#     Computational Statistics and Data Analysis, 44, 273-295
#
    
## Input:
##
## x0 - n x (p-1) matrix containing the explanatory variables;
## y  - n-vector containing binomial response (0 or 1);
    
## initwml  - logical value for selecting one of the two possible methods for computing
##          the initial value of the optimization process.
##          If initwml == TRUE (default), a weighted ML estimator is
##          computed with weights derived from the MCD estimator
##          computed on the explanatory variables. If initwml == FALSE,
##          a classical ML fit is perfomed.
##          When the explanatory variables contain binary observations,
##          it is recommended to set initwml to FALSE or to modify the
##          code of the algorithm to compute the weights only on the
##          continuous variables.
##
## const    - tuning constant used in the computation of the estimator (default=0.5);
## kmax     - maximum number of iterations before convergence (default=1000);
## maxhalf  - max number of step-halving (default=10).
##
## Output:
##
## A list with the follwoing components:
## convergence   - TRUE or FFALSE if convergence achieved or not
## objective     - value of the objective function at the minimum
## coef          - estimates for the parameters
## sterror       - standard errors of the parameters (if convergence is TRUE)
##
##Example:
##
## x0 <- matrix(rnorm(100,1))
## y  <- as.numeric(runif(100)>0.5)        #numeric(runif(100)>0.5)
## BYlogreg(x0,y)
##
    
    tmp=model.frame(formula, dat)
    y=tmp[,1]
    x0=model.matrix(formula, dat)
    if (colnames(x0)[1]=="(Intercept)") x0=x0[,-1,drop=FALSE] 
    
    y <- data.matrix(y)
    if(!is.numeric(y))
        y <- as.numeric(y)
    if(dim(y)[2] != 1)
        stop("y is not onedimensional")
    
        
    if(is.data.frame(x0))
    {
        x0 <- data.matrix(x0)
    }else if (!is.matrix(x0))
    {
        x0 <- matrix(x0, length(x0), 1,
                    dimnames = list(names(x0), deparse(substitute(x0))))
    }
    x <- as.matrix(cbind("Intercept"=1, x0))
    
    if(nrow(x) != nrow(y))
        stop("Number of observations in x and y not equal")
    
    na.x <- !is.finite(rowSums(x))
    na.y <- !is.finite(y)
    ok <- !(na.x | na.y)
    x <- x[ok, , drop = FALSE]
    y <- y[ok, , drop = FALSE]
    dx <- dim(x)
    n <- dx[1]
    if (n == 0)
        stop("All observations have missing values!")
    
    n <- nrow(x)
    p <- ncol(x)
    
    ## Smallest value of the scale parameter before implosion
    sigmamin <- 1e-4
    
    ## Computation of the initial value of the optimization process
#    if(initwml == TRUE)
#    {
#        hp <- floor(n*(1-0.25))+1
#    
###        mcdx <- cov.mcd(x0, quantile.used =hp,method="mcd")
###        rdx=sqrt(mahalanobis(x0,center=mcdx$center,cov=mcdx$cov))
#        mcdx <- rrcov::CovMcd(x0, alpha=0.75)
#        rdx  <- sqrt(rrcov::getDistance(mcdx))
#        vc   <- sqrt(qchisq(0.975, p-1))
#        wrd  <- rdx <= vc
#        gstart <- glm(y~x0, family=binomial, subset=wrd)$coef
#        names(gstart)[-1]=colnames(x0)
#    }else
#    {
        gstart <- glm(y~x0,family=binomial)$coef
        names(gstart)[-1]=colnames(x0)
#    }
    
    sigmastart=1/sqrt(sum(gstart^2))
    xistart=gstart*sigmastart
    stscores=x %*% xistart
    sigma1=sigmastart
    
    ## Initial value for the objective function
    oldobj <- mean(phiBY3(stscores/sigmastart,y,const))
    kstep <- jhalf <- 1
    
    while(kstep < kmax & jhalf < maxhalf)
    {
        unisig <- function(sigma)
        {
            mean(phiBY3(stscores/sigma,y,const))
        }
    
        optimsig=nlminb(sigma1,unisig,lower=0)
        sigma1=optimsig$par
    
        if(sigma1<sigmamin)
        {
            #print("Explosion")
            kstep=kmax
        }else
        {
            gamma1=xistart/sigma1
            scores=stscores/sigma1
            newobj=mean(phiBY3(scores,y,const))
            oldobj=newobj
            gradBY3=colMeans((derphiBY3(scores,y,const)%*%matrix(1,ncol=p))*x)
            h=-gradBY3+((gradBY3 %*% xistart) *xistart)
            finalstep=h/sqrt(sum(h^2))
            xi1=xistart+finalstep
            xi1=xi1/(sum(xi1^2))
            scores1=(x%*%xi1)/sigma1
            newobj=mean(phiBY3(scores1,y,const))
    
            ####stephalving
            hstep=jhalf=1
            while(jhalf <= maxhalf & newobj > oldobj)
            {
                hstep=hstep/2
                xi1=xistart+finalstep*hstep
                xi1=xi1/sqrt(sum(xi1^2))
                scores1=x%*%xi1/sigma1
                newobj=mean(phiBY3(scores1,y,const))
                jhalf=jhalf+1
            }
    
            if(jhalf == maxhalf+1 & newobj > oldobj)
            {
                if (verbose) print("Convergence Achieved")
            }else
            {
                jhalf=1
                xistart=xi1
                oldobj=newobj
                stscores=x%*% xi1
                kstep=kstep+1
            }
        }
    }
    
    gammaest=xistart/sigma1
    result <- list(convergence=FALSE, objective=0, coef=drop(gammaest), formula=formula)
    if(kstep < kmax) {
        stander=sterby3(x0, y, const, gammaest)
        result$convergence=TRUE
        result$objective=oldobj
        result$sterror=stander
    }    
    class(result)=c("rlogit",class(result))# don't add auc to this list, b/c it may be confusing to call auc.coef which will drop the intercept, 
    
    return(result)
}

## Functions needed for the computation of estimator of Bianco and Yohai
phiBY3 <- function(s,y,c3)
{
  s=as.double(s)
  dev=log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
  return(rhoBY3(dev,c3)+GBY3Fs(s,c3)+GBY3Fsm(s,c3))
}

rhoBY3 <- function(t,c3)
{
  (t*exp(-sqrt(c3))*as.numeric(t <= c3))+
    (((exp(-sqrt(c3))*(2+(2*sqrt(c3))+c3))-(2*exp(-sqrt(t))*(1+sqrt(t))))*as.numeric(t >c3))
}

psiBY3 <- function(t,c3)
{
    (exp(-sqrt(c3))*as.numeric(t <= c3))+(exp(-sqrt(t))*as.numeric(t >c3))
}

derpsiBY3 <- function(t,c3)
{
    res=NULL
    for(i in 1:length(t))
    {
        if(t[i] <= c3)
        {
            res=rbind(res,0)
        } else
        {
            res=rbind(res,-exp(-sqrt(t[i]))/(2*sqrt(t[i])))
        }
    }
    res
}

sigmaBY3 <- function(sigma,s,y,c3)
{
    mean(phiBY3(s/sigma,y,c3))
}

derphiBY3=function(s,y,c3)
{
    Fs= exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
    ds=Fs*(1-Fs)
    dev=log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
    Gprim1=log(1+exp(-abs(s)))+abs(s)*(s<0)
    Gprim2=log(1+exp(-abs(s)))+abs(s)*(s>0)

    return(-psiBY3(dev,c3)*(y-Fs)+((psiBY3(Gprim1,c3)-psiBY3(Gprim2,c3))*ds))
}

der2phiBY3=function(s, y, c3)
{
    s=as.double(s)
    Fs= exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
    ds=Fs*(1-Fs)
    dev=log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
    Gprim1=log(1+exp(-abs(s)))+abs(s)*(s<0)
    Gprim2=log(1+exp(-abs(s)))+abs(s)*(s>0)
    der2=(derpsiBY3(dev,c3)*(Fs-y)^2)+(ds*psiBY3(dev,c3))
    der2=der2+(ds*(1-2*Fs)*(psiBY3(Gprim1,c3)-psiBY3(Gprim2,c3)))
    der2=der2-(ds*((derpsiBY3(Gprim1,c3)*(1-Fs))+(derpsiBY3(Gprim2,c3)*Fs)))

    der2
}


GBY3Fs <- function(s,c3)
{
    Fs= exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
    resGinf=exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fs))))-1)
    resGinf=(resGinf+(Fs*exp(-sqrt(-log(Fs)))))*as.numeric(s <= -log(exp(c3)-1))
    resGsup=((Fs*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s > -log(exp(c3)-1))

    return(resGinf+resGsup)
}


GBY3Fsm <- function(s,c3)
{
  Fsm=exp(-(log(1+exp(-abs(s)))+abs(s)*(s>0)))
  resGinf=exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fsm))))-1)
  resGinf=(resGinf+(Fsm*exp(-sqrt(-log(Fsm)))))*as.numeric(s >= log(exp(c3)-1))
  resGsup=((Fsm*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s < log(exp(c3)-1))
  return(resGinf+resGsup)
}

## Compute the standard erros of the estimates -
##  this is done by estimating the asymptotic variance of the normal
##  limiting distribution of the BY estimator - as derived in Bianco
##  and Yohai (1996)
##
sterby3 <- function(x0, y, const, estim)
{
    n <- dim(x0)[[1]]
    p <- dim(x0)[[2]]+1

    z <- cbind(rep(1,n),x0)
    argum <- z %*% estim

    matM <- matrix(data=0, nrow=p, ncol=p)
    IFsquar <- matrix(data=0, nrow=p, ncol=p)
    for(i in 1:n)
    {
        myscalar <- as.numeric(der2phiBY3(argum[i],y[i],const))
        matM <- matM + myscalar * (z[i,] %*% t(z[i,]))
        IFsquar <- IFsquar + derphiBY3(argum[i],y[i],const)^2 * (z[i,] %*% t(z[i,]))
    }

    matM    <- matM/n
    matMinv <- solve(matM)
    IFsquar <- IFsquar/n
    asvBY   <- matMinv %*% IFsquar %*% t(matMinv)

    sqrt(diag(asvBY))/sqrt(n)
}



coef.rlogit=function(object, ...) object$coef

ratio.rlogit=ratio.auc

predict.rlogit=function(object, newdata, ...) {
    predict.auc (object, newdata)        
}

trainauc.rlogit=function(fit, training.data=NULL, ...) {
    trainauc.auc(fit, training.data)
}
