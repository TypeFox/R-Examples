"pois.exact" <-
function(x, pt = 1, conf.level = 0.95){
  ## x = Poisson count
  ## pt = person time
  ##updated 2004-11-29
  xc <- cbind(x,conf.level,pt)
  pt2 <- xc[,3]
  results <- matrix(NA,nrow(xc),6)
  f1 <- function(x,ans,alpha=alp) {ppois(x,ans)-alpha/2}
  f2 <- function(x,ans,alpha=alp) 1-ppois(x,ans)+dpois(x,ans)-alpha/2
  for(i in 1:nrow(xc)){
    alp <- 1-xc[i,2]
    interval <- c(0,xc[i,1]*5+4)
    uci <- uniroot(f1,interval=interval,x=xc[i,1])$root/pt2[i]
    if(xc[i,1]==0){
      lci <- 0
    } else {
      lci <- uniroot(f2,interval=interval,x=xc[i,1])$root/pt2[i]
    } 
    results[i,] <- c(xc[i,1],pt2[i],xc[i,1]/pt2[i],lci,uci,xc[i,2]) 
  }
  coln <- c("x","pt","rate","lower","upper","conf.level")
  colnames(results) <- coln
  data.frame(results)
}
