SubsCV <- function(x,MaxLag,parallel){
 n <- length(x)
 if (((n-MaxLag)<0) || ((n-MaxLag)<4) || ((n-MaxLag)<=25)) stop("Give bigger sample size n")
 test <- function(x,j,b){
  q <- n-b+1 
  q <- ifelse(q<0,abs(q),q)
  Xb <- sapply(1:q, function(i) x[i:(i+b-1)])
  Vx <- unlist(lapply(1:q, FUN=function(i) ADCF(Xb[,i],MaxLag=j)[j+1]^2))
  Vstar <- sort(Vx)
  la <- ceiling((1-0.05)*q)
  critical.value <- sqrt((b-j)*Vstar[la]/(n-j))
  return(critical.value)
 }
 optimal.block <- function(x,lags,MaxLag){
  r <- 2
  bsmall <- MaxLag+4
  bbig <- bsmall+20
  bseq <- seq(bsmall,bbig,by=4)
  se <- rep(0,length(bseq))
  for (l in 1:length(bseq)){
   bk <- bseq[l]+seq(-r,r,by=1)
   Tneighb <- sapply(1:length(bk), function(m) test(x,j=lags,b=bk[m]))  
   se[l] <- (sum(abs(Tneighb-mean(Tneighb))^2)/(length(bk)-1))^(1/2)
  }
  smallest <- order(se)[1]
  b.star <- bseq[smallest]
  return(b.star)
 }
 if(parallel==TRUE){
  closeAllConnections()
  #cl <- makeCluster(detectCores())
  cl <- makeCluster(2)
  registerDoParallel(cl)
  clusterSetRNGStream(cl = cl, iseed = 9182)
  i <- 1:MaxLag
  fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="c",.export=c("ADCF","dcor")) ))
  fe <- eval(fe_call)
  cv <- fe %dopar% test(x,i,b=optimal.block(x,lags=i,MaxLag))
  stopCluster(cl)
 } else {
 cv=sapply(1:MaxLag, function(i) test(x,i,b=optimal.block(x,lags=i,MaxLag)))
 }
 return(cv)
}

