# simulating data from a DCCC-GARCH(1,1) and computing a 
# DCCC-GARCH(1,1) volatility and dynamic conditional correlations


dcc.sim <- function(nobs, a, A, B, R, dcc.para, d.f = Inf, cut=1000, model){
  if(model=="diagonal"){
      A <- diag(diag(A)); B <- diag(diag(B))
  }
  if(stationarity(A, B)>1){
      stop("A and B matrices do not satisfy the stationarity condition \n")
  }
  nobs <- nobs+cut
  Id <- diag(length(a))
  inih <- solve(Id-A-B)%*%a
  dccpar1 <- dcc.para[1]; dccpar2 <- dcc.para[2]
  sim <- .Call("dcc_sim", nobs, a, A, B, inih, R, dccpar1, dccpar2, d.f)
  list(  z     = sim[[1]][(cut+1):(nobs),],
         std.z = sim[[2]][(cut+1):(nobs),],
         dcc   = sim[[3]][(cut+1):(nobs),],
         h     = sim[[4]][(cut+1):(nobs),],
         eps   = sim[[5]][(cut+1):(nobs),])
}
