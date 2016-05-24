ci.p<-function(data,conf=.95,summarized=FALSE, phat=NULL,S.phat=NULL,fpc=FALSE,n=NULL,N=NULL,method="agresti.coull",plot=TRUE){
indices <- c("agresti.coull","asymptotic","exact","LR","score")
method <- match.arg(method, indices)
alpha <-1-conf 
z.star<-qnorm(1-(alpha/2)) 

if(summarized==FALSE){
  n<-length(as.matrix(data))
  phat<-ifelse(method == "agresti.coull", (sum(data) + 2)/(n + 4), sum(data)/n)
  Var.phat<-ifelse(method == "agresti.coull", (phat * (1 - phat))/(n + 4), (phat * (1 - phat))/n)
  S.phat<-ifelse(fpc == FALSE, sqrt(Var.phat), sqrt((1-(n/N)) * Var.phat))
}

if(summarized==TRUE){
  S.phat<-ifelse(fpc==FALSE,S.phat,sqrt((1-(n/N))*S.phat))}
  
  m <- S.phat * z.star
    
if(method=="agresti.coull"|method=="asymptotic")CI<-c(phat,phat-m,phat+m)
  
  x <- n * phat

if(method == "exact") {
   cl <- qbeta(alpha/2, x, n - x + 1)
	 cu <- qbeta(1 - (alpha/2), x + 1, n - x)	
	 CI <- c(phat, cl, cu)
   }

if(method=="score"){ 
  a <- 1 + (z.star^2/n)
  b <- -2*phat - (z.star^2/n)
  c <- phat^2 
  CI<-c(phat,(-b-sqrt(b^2-(4*a*c)))/(2*a),(-b+sqrt(b^2-(4*a*c)))/(2*a)) # quadratic eqn. 
  }

if(method=="LR"){ 
  x2 <- qchisq(conf, 1)/2
  p.0 <- seq(0.00001, 0.99999, by = 0.00001)
  LR <- rep(NA, length=length(p.0))
  ML <- log(dbinom(x,n,phat))
  LL <- log(dbinom(x,n,p.0))
  int<-p.0[LL > ML - x2]
  CI <- c(phat, min(int), max(int))
  if(plot == TRUE){
  plot(p.0, LL, type="l", xlab = expression(pi[0]), ylab = "Binomial log-likelihood") 
  abline(h = ML-x2, lty=2); abline(v=CI,lwd=2, col="gray"); abline(v=phat, lwd = 2)
  #legend("bottomleft",lty=2,legend=expression(paste("logLik(",hat(pi),") - ",italic(f),"(",pi[0],") = ",chi^2,"(1 - ",alpha,")/2", sep = "")), box.col = "white", cex = .9, bg = "white", inset = .01)
  } 
  } 
  

head<-paste(paste(as.character(conf*100),"%",sep=""),c("Confidence interval for binomial parameter pi"))
if(method=="agresti.coull")head<-paste(head,"(method=Agresti-Coull)")
if(method=="exact")head<-paste(head,"(method=Clopper-Pearson)")
if(method=="asymptotic")head<-paste(head,"(method= asymptotic normal approximation)")
if(method=="score")head<-paste(head,"(method=score)")
if(method=="LR")head<-paste(head,"(method=likelihood ratio)")


ends<-c("Estimate",paste(as.character(c((1-conf)/2,1-((1-conf)/2))*100),"%",sep=""))
if(method=="agresti.coull"|method=="asymptotic"|method =="score"|method=="LR")res<-list(SE=S.phat,margin=m,ci=CI,ends=ends,head=head)
if(method=="exact")res<-list(ci=CI,ends=ends,head=head)
class(res)<-"ci"
res
}
