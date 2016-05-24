dw.reg<-function(formula, data,tau=0.5,para.q1=FALSE,para.q2=TRUE,para.beta=FALSE, ...)
{
loglik.q1<-function(par,x,y)
{
beta<-par[length(par)]
if(beta<=0)
	loglik<-NA
else
{
theta<-par[-length(par)]
logq<-x%*%theta-log(1+exp(x%*%theta))
if(any(exp(logq)==0) | any(exp(logq)==1))
	loglik<-NA
else
	loglik<-sum(log(exp(y^beta*logq)-exp((y+1)^beta*logq)))
}
return(loglik)
}

loglik.q2<-function(par,x,y)
{
  beta<-par[length(par)]
  if(beta<=0)
    loglik<-NA
  else
  {
    theta<-par[-length(par)]
    logq<- -exp(x%*%theta)
    if(any(exp(logq)==0) | any(exp(logq)==1))
      loglik<-NA
    else
      loglik<-sum(log(exp(y^beta*logq)-exp((y+1)^beta*logq)))
  }
  return(loglik)
}


loglik.beta<-function(par,x,y)
{
q<-par[length(par)]
if(q<=0 | q >=1)
	loglik<-NA
else
{
theta<-par[-length(par)]
beta<-exp(x%*%theta)
loglik<-sum(log(exp(y^beta*log(q))-exp((y+1)^beta*log(q))))
}
return(loglik)
}


loglik.q1beta<-function(par,x,y)
{
theta.q<-par[1:ncol(x)]
theta.beta<-par[(ncol(x)+1) : length(par)]
beta<-exp(x%*%theta.beta)
logq<-x%*%theta.q-log(1+exp(x%*%theta.q))
if(any(exp(logq)==0) | any(exp(logq)==1))
	loglik<-NA
else
	loglik<-sum(log(exp(y^beta*logq)-exp((y+1)^beta*logq)))
return(loglik)
}

loglik.q2beta<-function(par,x,y)
{
  theta.q<-par[1:ncol(x)]
  theta.beta<-par[(ncol(x)+1) : length(par)]
  beta<-exp(x%*%theta.beta)
  logq<--exp(x%*%theta.q)
  if(any(exp(logq)==0) | any(exp(logq)==1))
    loglik<-NA
  else
    loglik<-sum(log(exp(y^beta*logq)-exp((y+1)^beta*logq)))
  return(loglik)
}


call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
y <- model.response(mf, "numeric")
x <- model.matrix(mt, mf, contrasts)
term.labels <- colnames(x)

oldw <- getOption("warn")
options(warn = -1)

q<-try(estdweibull(y,method="ML",zero=TRUE),silent=TRUE)[1]
beta<-try(estdweibull(y,method="ML",zero=TRUE),silent=TRUE)[2]

if(is.na(q) & para.beta)
{
  q<-try(estdweibull(y,method="P",zero=TRUE),silent=TRUE)[1]
}
if(is.na(q) & para.beta)
{
  q<-try(estdweibull(y,method="M",zero=TRUE),silent=TRUE)[1]
}
if(is.na(q) & para.beta)
{
  set.seed(131)
  q<-runif(1)
}

if(is.na(beta) & (para.q1 | para.q2))
{
  beta<-try(estdweibull(y,method="P",zero=TRUE),silent=TRUE)[2]
}
if(is.na(beta) & (para.q1 | para.q2))
{
  beta<-try(estdweibull(y,method="M",zero=TRUE),silent=TRUE)[2]
}
if(is.na(beta) & (para.q1 | para.q2))
{
  set.seed(131)
  beta<-runif(1)
}
options(warn = oldw)

if(para.q1){
int.in<-log(q/(1-q))
par.in<-c(int.in,rep(0,ncol(x)-1),beta)
names(par.in)<-c(term.labels,"beta")
mle<-maxLik(loglik.q1,start=par.in,x=x,y=y,...)
est<-coef(mle)
beta<-est[length(est)]
theta<-est[-length(est)]
logq<-x%*%theta-log(1+exp(x%*%theta))
Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/beta)-1)
}

if(para.q2){
  int.reg<--1*as.numeric(coef(glm(y~x-1,family=quasipoisson)))
  par.in<-c(int.reg,beta)
  names(par.in)<-c(term.labels,"beta")
  mle<-maxLik(loglik.q2,start=par.in,x=x,y=y,...)
  est<-coef(mle)
  beta<-est[length(est)]
  theta<-est[-length(est)]
  logq<--exp(x%*%theta)
  Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/beta)-1)
}


if(para.beta)
{
int.in<-log(beta)
par.in<-c(int.in,rep(0,ncol(x)-1),q)
names(par.in)<-c(term.labels,"q")
mle<-maxLik(loglik.beta,start=par.in,x=x,y=y,...)
est<-coef(mle)
q<-est[length(est)]
theta<-est[-length(est)]
logbeta<-x%*%theta
Fittedy.tau<-ceiling((log(1-tau)/log(q))^(1/exp(logbeta))-1)
}

if(para.beta & para.q1)
{
int.in.beta<-log(beta)
int.in.q<-log(q/(1-q))
par.in<-c(int.in.q,rep(0,ncol(x)-1),int.in.beta,rep(0,ncol(x)-1))
names(par.in)<-c(term.labels,term.labels)
mle<-maxLik(loglik.q1beta,start=par.in,x=x,y=y,...)
est<-coef(mle)
theta.q<-est[1: (ncol(x))]
theta.beta<-est[(ncol(x)+1) : length(est)]
logq<-x%*%theta.q-log(1+exp(x%*%theta.q))
logbeta<-x%*%theta.beta
Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/exp(logbeta))-1)
}


if(para.beta & para.q2)
{
  int.in.beta<-log(beta)
  int.in.q<-log(-log(q))
  par.in<-c(int.in.q,rep(0,ncol(x)-1),int.in.beta,rep(0,ncol(x)-1))
  names(par.in)<-c(term.labels,term.labels)
  mle<-maxLik(loglik.q2beta,start=par.in,x=x,y=y,...)
  est<-coef(mle)
  theta.q<-est[1: (ncol(x))]
  theta.beta<-est[(ncol(x)+1) : length(est)]
  logq<--exp(x%*%theta.q)
  logbeta<-x%*%theta.beta
  Fittedy.tau<-ceiling((log(1-tau)/logq)^(1/exp(logbeta))-1)
}
fit <- list()
fit$call <- call
fit$coefficients <- coef(mle)
fit$loglik <- logLik(mle)
fit$fitted.values <- Fittedy.tau
fit$tTable <- summary(mle)
class(fit) <- "dw.reg"
return(fit)
}
