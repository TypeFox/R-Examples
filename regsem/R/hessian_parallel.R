
hessian_parallel = function(par,ImpCov,A,S,F,A_fixed,A_est,S_fixed,S_est,lambda,alpha,type){

  hesS_out <- matrix(0,length(par),length(par))
  h <- 0.0001

  # using Cudeck, Klebe, Henly (1993)

  add <- matrix(0,length(par),length(par))
  diag(add) <- h


  vec = seq(1:length(par))
  li <- list()
  grid <- expand.grid(i=vec,j=vec)
  grid <- as.matrix(grid)
  for(i in 1:nrow(grid)) li[[i]] <- grid[i,]

  snowfall::sfExport("add","par","A","ImpCov","A_fixed","A_est","S","S_fixed","S_est","F","li")


  hess_fun <- function(indexI,indexJ,ImpCov,A,S,F,A_fixed,A_est,S_fixed,S_est,add){

    ImpCovI = RAMmult((par + add[indexI,]),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
    ImpCovJ = RAMmult((par + add[,indexJ]),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
    ImpCovII <- (ImpCovI - ImpCov)/h
    ImpCovJJ <- (ImpCovJ - ImpCov)/h

     0.5 * trace(solve(ImpCov) %*% ImpCovII %*% solve(ImpCov) %*% ImpCovJJ)
  }




  matrix(unlist(snowfall::sfLapply(li,function(x) hess_fun(indexI=x["i"],indexJ=x["j"],
                ImpCov,A,S,F,A_fixed,A_est,S_fixed,S_est,add))),length(par),length(par))


}

