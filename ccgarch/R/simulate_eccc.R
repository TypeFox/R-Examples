# simulating data from an ECCC-GARCH(1,1)


   eccc.sim <- function(nobs, a, A, B, R, d.f=Inf, cut=1000, model){
     if(model=="diagonal"){
         A <- diag(diag(A)); B <- diag(diag(B))
     }
     if(stationarity(A, B)>1){
         stop("A and B matrices do not satisfy the stationarity condition. \n")
     }
      nobs <- nobs+cut
      Id <- diag(length(a))
      inih <- solve(Id-A-B)%*%a
#      sim <- .Call("eccc_sim", nobs, a, A, B, R, inih, d.f, PACKAGE="ccgarch")
      sim <- .Call("eccc_sim", nobs, a, A, B, R, inih, d.f)
      list( h = sim[[1]][(cut+1):(nobs),],
            eps = sim[[2]][(cut+1):(nobs),])
   }
