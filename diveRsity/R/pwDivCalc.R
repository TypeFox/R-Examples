################################################################################
# pwDivCalc: a small function for calculating pairwise ht and hs 
################################################################################
pwDivCalc <- function(x, pw, npops){
  ht <- matrix(ncol = npops, nrow = npops)
  hs <- matrix(ncol = npops, nrow = npops)
  for(i in 1:ncol(pw)){
    gamma <- sum(sqrt(abs(x[,pw[1,i]] * x[,pw[2,i]])))^-1 
    f <- gamma * sqrt(x[,pw[1,i]] * x[,pw[2,i]])
    ht[pw[1,i],pw[2,i]] <- 1 - sum(((f + x[,pw[1,i]])/2)^2)
    ht[pw[2,i],pw[1,i]] <- 1 - sum(((f + x[,pw[2,i]])/2)^2)
    hs[pw[1,i],pw[2,i]] <- 1 - sum((f^2 + x[,pw[1,i]]^2)/2)
    hs[pw[2,i],pw[1,i]] <- 1 - sum((f^2 + x[,pw[2,i]]^2)/2)
  }
  ht[is.nan(ht)] <- 0
  hs[is.nan(hs)] <- 0
  list(ht = ht, 
       hs = hs)
}
################################################################################
# END - pwDivCalc
################################################################################