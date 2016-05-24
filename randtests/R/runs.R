##
##  probability function of the runs statistic
##
druns <- function(x, n1, n2, log = FALSE){
  stopifnot(is.numeric(x))
  x <- ifelse(x == round(x),x,1)
  r0 <- ifelse(x %% 2==0, 2*choose(n1-1, round(x/2)-1)*choose(n2-1, round(x/2)-1), 
             choose(n1-1, round((x-1)/2))*choose(n2-1, round((x-3)/2))+choose(n1-1, round((x-3)/2))*choose(n2-1, round((x-1)/2)))  
  r0<-r0/choose(n1+n2, n1)
# if TRUE, probabilities p are given as log(p).  
ifelse(log,return(log(r0)),return(r0))  
}

##
##  distribution function of the runs statistic
##
pruns <- function(q, n1, n2, lower.tail = TRUE, log.p = FALSE){
  stopifnot(is.numeric(q) & n1>0 & n2>0)
  q <- ifelse(q >= 1, q, 1)
  q <- ifelse(q <= n1+n2, q, n1+n2)
  q <- round(q)
  tmp <- cumsum(druns(1:max(q),n1,n2,log=log.p))
  r0 <- tmp[q]
  if (lower.tail==FALSE){r0<- 1-r0}  
#  r0 <- NULL
#  if (lower.tail){
#    for (i in 1:length(q)){r0 <- c(r0,ifelse(q[i]>=2,sum(druns(x=2:floor(q[i]),n1,n2,log=log)),0))}
#  }
#  else {r0 <- 1-pruns(q,n1,n2,lower.tail=T, log=log)}  
  return(r0)  
}  
##
##  Quantile function of the runs statistic
##
qruns <- function(p, n1, n2, lower.tail = TRUE, log.p = FALSE){
  r0 <- NULL
  q1 <- ifelse (n1==n2, 2*n1, 2*min(n1,n2)+1) 
  pr <- c(0, cumsum(druns(2:q1, n1, n2)))
  for (i in 1:length(p)){
    if (p[i]>=0 & p[i]<=1){
      #rq<-which(abs(pr-p)==min(abs(pr-p))) 
      qr <- NULL
      for (j in 2:q1){
        if (pr[j-1]<p[i] & p[i]<=pr[j]){qr<-j}
      }
      if (p[i] == pr[1]){qr <- 2}
    }
    else {rq<-NA}
    r0<-c(r0, qr)
  }
  return(r0)  
}  
##
##  Generates (pseudo) randon values of the runs statistic
##
rruns <- function(n, n1, n2){
  return(qruns(runif(n), n1, n2))  
}    
