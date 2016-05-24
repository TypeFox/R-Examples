
WNS.test <- function(x,type="TS",s=0,r=0, nsim=1000){
  MLE<-function(x){
    if(sum(x<0)){stop(paste("Data x is not a positive sample"))}
    n=length(x)
    #Apply the log. transformation => If: X -> Weibull, then -log(X) -> EV
    lv = -log(x)
    y = sort(lv)
    #Compute the MLEs of the EV distribution
    f1 <- function(toto,vect){
      if(toto!=0){
        f1= sum(vect)/length(vect)
        f1 = f1 - sum(vect*exp(-vect*toto))/sum(exp(-vect*toto)) - 1/toto
      }else{ f1=100}
      f1=abs(f1)
    }  
    theta <- optimize(f1,c(0.0001,50),maximum=FALSE,vect=y,tol=1/10^5)
    t <- theta$minimum
    aux <- sum(exp(-y*t)/length(x))
    ksi <- -(1/t)*log(aux)
    #Compute the pseudo-observation y_1, .., y_n
    y <- -(y-ksi)*t
    return(MLE<-list(eta=exp(-ksi),beta=t,y=y))
  }
  ##Family of the test statistics based on the normalized spacings
  NS.statistic<-function(x,type,s,r,X=NULL){
    TYPE <- deparse(substitute(type))
    m <- length(x)
    n <- m+s+r
    A <- sort(log(x))
    d1 <- A[2:(m-1)]-A[1:(m-2)]
    d2 <- A[2:m]-A[1:(m-1)]
    if(is.null(X)){X <- GoFNS(r+1,n,m)}
    mu1 <- X[2:(m-1)]-X[1:(m-2)]
    mu2 <- X[2:m]-X[1:(m-1)]
    l <- seq(length=(m-2),from=(r+1),to=(n-s-2))
    G1 <- d1/mu1
    G2 <- d2/mu2
    w1 <- 2*(sum((n-s-1-l)*G1))
    w2 <- (m-2)*sum(G2)
    #Tiku-Singh test statistic 
    if(type=="TS"){
      NS.statistic <- w1/w2
      #Lockhart-O'Reilly-Stephens test statistic 
    }else if(type=="LOS"){
      z <- NULL
      for (i in 1:(m-2)){z[i]=sum(G2[1:i])/sum(G2)}
      z <- sort(z)
      z1 <- sort(z,decreasing=TRUE)
      I <- 1:(m-2)
      NS.statistic <- -(m-2)-1/(m-2)*sum((2*I-1)*(log(z)+log(1-z1)))
      #Mann-Scheuer-Fertig statistic only right censoring
    }else if(type=="MSF"){
      if(s!=0){stop(paste('the test is only applied for right censoring'))}
      l1 <- floor(m/2)
      l2 <- m-l1-1
      S <- sum( (A[(l1+2):m]-A[(l1+1):(m-1)]) / (X[(l1+2):m]-X[(l1+1):(m-1)]) )
      S <- S/sum((A[2:m]-A[1:(m-1)])/(X[2:m]-X[1:(m-1)]))
      #  L <- l2/l1*sum((A[(l1+2):m]-A[(l1+1):(m-1)])/(X[(l1+2):m]-X[(l1+1):(m-1)]))
      #  L <- L/sum((A[2:(l2+1)]-A[1:(l2)])/(X[2:(l2+1)]-X[1:l2])) 
      NS.statistic <- S#,L)
    }
    return(NS.statistic)
  }
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="TS"){
    METHOD="Test of Tiku and Singh for the Weibull distribution"
  } else if(as.character(type)=="LOS"){
    METHOD="Test of Lockhart-O'Reilly and Stephens for the Weibull distribution"
  } else if(as.character(type)=="MSF"){
    METHOD="Test of Mann-Scheuer and Fertig for the Weibull distribution"}
  MLE <- MLE(x)
  estimate.obs <- c(MLE$eta,MLE$beta)
  m <- length(x)
  n <- m+s+r
  X <- GoFNS(r+1,n,m)
  statistic.obs <- NS.statistic(x,type,s,r,X)  
  fun<-function(y){
    y<-sort(y)[(r+1):(n-s)]
    fun <- NS.statistic(y,type,s,r,X)
     
    return(fun)
  }
  sim.statistic <- GoFsim(nsim,n,fun)
  sim.statistic <- sort(sim.statistic)
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  pvalb <- 2*min(p_val,1-p_val)
  p_val <- switch(type,"TS"=pvalb,"MSF"=pvalb,p_val)
  WNS.test <- list(statistic =c(S=statistic.obs), p.value = p_val, 
                   method=paste(METHOD),
                   estimate = c(eta=estimate.obs[1],beta=estimate.obs[2]),  
                   data.name=DNAME)
  
  class(WNS.test) <-"htest"
  return(WNS.test)
}
