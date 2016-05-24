dsample.pps <-
function(mos,n)
{ 
  ## Description: sample size with probability
  ##		proportional to measure of size
  ## mos : measure of sizes each unit
  ##  n  : sample size
  ## value: sample sizes for each unit
  if(all(mos==0) | n==0) return(rep(0,length(mos)))
  nraw <- n * mos / sum(mos)
  nn <- length(mos)
  n1 <- floor(nraw)
  p <- nraw - n1
  nfrac <- sum(p)
  if ( nfrac >0){
  	addoneid <- sample(nn,nfrac,prob=p)
  	n1[addoneid] <- n1[addoneid] + 1
  }
  ##stopifnot(sum(n1)==n)
  n1
}
