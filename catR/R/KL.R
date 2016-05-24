KL <- function (itemBank, item, x, it.given, model=NULL, theta=NULL,lower = -4, upper = 4, nqp = 33, 
                type = "KL", priorDist="norm", priorPar = c(0, 1), D=1, X=NULL, lik = NULL) 
{
  if (type != "KL" & type != "KLP") 
    stop("'type' must be either 'KL' or 'KLP'", call. = FALSE)
if (!is.null(X) & !is.null(lik)){
if (length(X) != length(lik)) stop("'X' and 'lik' must have the same length!",call.=FALSE)
}
if (is.null(theta)) theta<-thetaEst(it.given,x,D=D,model=model, method="ML")
  KLF<-NULL
  par <- rbind(itemBank[item,])
if (is.null(X)) X<-seq(from=lower,to=upper,length=nqp)
if (is.null(model)){
if (is.null(lik)){
L <- function(th, r, param) 
      prod(Pi(th, param,D=D)$Pi^r * (1 - Pi(th,param,D=D)$Pi)^(1 - r))
lik<-sapply(X,L,x,it.given)
}
  KLF[1:nqp] <- Pi(theta,par,D=D)$Pi * log(Pi(theta,par,D=D)$Pi/Pi(X[1:nqp],par,D=D)$Pi) + (1 - Pi(theta,par,D=D)$Pi) * log((1 - Pi(theta,par,D=D)$Pi)/(1 - Pi(X[1:nqp],par,D=D)$Pi))
  crit.value <- lik*KLF 
  if (type=="KLP") {
    pd<-switch(priorDist,norm=dnorm(X,priorPar[1],priorPar[2]),unif=dunif(X,priorPar[1],priorPar[2]))
    crit.value <- crit.value*pd 
  }
}
else{
if (is.null(lik)){
LL <- function(th, param, r, model,D=1) {
prob <- Pi(th, param, model = model,D=D)$Pi
res <- 1
for (i in 1:length(r)) res <- res * prob[i, r[i] + 1]
return(res)
}
lik<-sapply(X,LL,it.given,x,model=model,D=D)
}
pi<-Pi(theta,par,model=model,D=D)$Pi
for (i in 1:length(X)){
pri<-Pi(X[i],par,model=model,D=D)$Pi
KLF[i]<-sum(pi*log(pi/pri),na.rm=TRUE)
}
crit.value <- lik*KLF 
if (type=="KLP") {
pd<-switch(priorDist,norm=dnorm(X,priorPar[1],priorPar[2]),unif=dunif(X,priorPar[1],priorPar[2]))
crit.value <- crit.value*pd 
}
}
RES <- integrate.catR(X, crit.value)
return(RES)
}
