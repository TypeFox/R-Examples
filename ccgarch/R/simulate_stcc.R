# simulating data from a STCC-GARCH(1,1) and computing a 
# STCC-GARCH(1,1) volatility and dynamic conditional correlations

#   source("st_functions.r")
   
   stcc.sim <- function(nobs, a, A, B, R1, R2, tr.par, st.par, d.f = Inf, cut=1000, model){
        if(model=="diagonal"){
            A <- diag(diag(A)); B <- diag(diag(B))
        }
        if(stationarity(A, B)>1){
            stop("A and B matrices do not satisfy the stationarity condition \n")
        }
      nobs <- nobs+cut; ndim <- length(a)
      tr.var <- uni.vola.sim(tr.par, nobs, d.f = Inf, cut = cut)$eps
      st <- tr.func(st.par, tr.var)
      vecR <- (1-st)*matrix(R1, nobs, ndim^2, byrow=TRUE) + st*matrix(R2, nobs, ndim^2, byrow=TRUE)
      Id <- diag(length(a))
      inih <- solve(Id-A-B)%*%a
#      sim <- .Call("stcc_sim", nobs, a, A, B, vecR, inih, d.f, PACKAGE="ccgarch")
      sim <- .Call("stcc_sim", nobs, a, A, B, vecR, inih, d.f)
      list( h = sim[[1]][(cut+1):nobs,],
            eps = sim[[2]][(cut+1):nobs,],
            tr.var = tr.var,
            st = st, vecR=vecR)
   }
