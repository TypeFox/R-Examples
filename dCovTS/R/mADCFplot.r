mADCFplot <- function(x,MaxLag=15,ylim=NULL,b=499){
if(b<=0) stop("No plots are given for b<=0")
if (MaxLag==0) stop("MaxLag must be greater than 0")
q <- as.integer(NCOL(x))
n <- as.integer(NROW(x))
R <- array(unlist(lapply(0:MaxLag, FUN=function(i) mADCF(x,i,output=FALSE))),dim=c(q,q,MaxLag+1))
mcv <- max(mRbootCV(n,q,MaxLag,b=b,parallel=TRUE))
Lag <- 0:MaxLag
Rplot <- function(R,cv,main){
 r1 <- max(max(cv),1)
 if (is.null(ylim)) ylim=c(0,r1)
 plot(Lag,R,type="n",main=main,ylab="ADCF",ylim=ylim)
 for (i in seq(0,MaxLag,by=1)){
  segments(i,0,i,R[(i+1)])
 }
 points(0:MaxLag,rep(cv,MaxLag+1),type="l",lty=3,lwd=2,col="blue")
}
plot.fun <- function(d1,d2){
plot.new()
for (i in 1:length(d1)){
 for (j in 1:length(d2)){
   par(mfg=c(i,j))
    res <- sapply(0:MaxLag,FUN=function(k) R[,,(k+1)][d1[i],d2[j]])
    if(length(z)==0){
      if (d1[i]==d2[j]) Rplot(res,mcv,paste("Series",d1[i]))
      else Rplot(res,mcv,paste("Series",d1[i],"&",d2[j]))
    } else {
      if (d1[i]==d2[j]) Rplot(res,mcv,paste(z[d1[i]]))
      else Rplot(res,mcv,paste(z[d1[i]],"&",z[d2[j]]))
    }
 }
}
}
if(q>8) {
 print("No plots are given due to high dimension")
}
else {
d <- ifelse(q>4, 4, q)
z <- colnames(x)
counter <- 1
s1 <- seq(1,d,by=1)
 if(q>4) {
   s2 <- seq(5,q,by=1)
 } else {
   s2 <-  NA
 }
 if (q<=4){
  par(mfrow=c(q,q))
  plot.fun(s1,s1)
 } else {
  for (sc in 1:4){
        #par(ask=TRUE)
   if (counter==1){
    par(mfrow=c(4,4))
    plot.fun(s1,s1)
   } else if (counter==2){
    par(ask=TRUE)
    par(mfrow=c(4,4))
    plot.fun(s1,s2)
   } else if (counter==3){
     par(mfrow=c(4,4))
     plot.fun(s2,s1)
    } else if (counter==4){
      par(mfrow=c(4,4))
      plot.fun(s2,s2)
     }
  counter <- counter+1
  }
}
}
result <- list(matrices=R,critical.value=mcv)
return(result)
}

