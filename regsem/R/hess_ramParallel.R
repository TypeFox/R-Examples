
hess_ramParallel = function(par,ImpCov,SampCov,Areg,Sreg,A,S,F){

  #?sfClusterApplyLB
  #sfClusterApplyLB(index, simulationIteration, n, rho, g, h, pv, het=F)

  hess.out <- matrix(0,length(par),length(par))



  # get rid of
  #Areg = A.est
  #Sreg= S.est
  #ImpCov = out$imp_cov
  #SampCov=S1



  B = solve(diag(nrow(A)) - Areg)
  C = diag(nrow(ImpCov)) - solve(ImpCov) %*% SampCov
  E = B %*% Sreg %*% t(B)


  vec = seq(1:length(par))
  li <- list()
  grid <- expand.grid(i=vec,j=vec)
  grid <- as.matrix(grid)
  for(i in 1:nrow(grid)) li[[i]] <- grid[i,]


  snowfall::sfExport("Areg","A","ImpCov","SampCov","B","C","E","S","li")

  hess_fun <- function(indexI,indexJ,Areg,A,ImpCov,SampCov,B,C,E,S){
    #indexI <- list[1]
    #indexJ <- list[2]
    Ai <- A == indexI
    Aj <- A == indexJ;
    Aij <- matrix(0,nrow(A),ncol(A))
    Aij[Ai==T | A==T] <- 1

    Si <- S == indexI;
    Sj <- S == indexJ;
    Sij <- matrix(0,nrow(S),ncol(S))
    Sij[Si==T | Sj==T] <- 1

    deriv15_I <- F %*% B %*% Ai %*% E %*% t(F) + F %*% B %*% Si %*% t(B) %*% t(F)
    deriv15_J <- F %*% B %*% Aj %*% E %*% t(F) + F %*% B %*% Sj %*% t(B) %*% t(F)

    deriv15_IJ <- F %*% B %*% Ai %*% B %*% Aj %*% E %*% t(F) +
      F %*% B %*% Aj %*% B %*% Ai %*% E %*% t(F) +
      F %*% B %*% Ai %*% B %*% Sj %*% t(B) %*% t(F) +
      F %*% B %*% Ai %*% E %*% t(Aj) %*% t(B) %*% t(F) +
      F %*% B %*% Aj %*% B %*% Si %*% t(B) %*% t(F)

    # left out mean part
   trace(solve(ImpCov) %*% deriv15_IJ %*% C - solve(ImpCov) %*%
                              deriv15_I %*% solve(ImpCov) %*% deriv15_J %*% C +
                              solve(ImpCov) %*% deriv15_J %*% solve(ImpCov %*% deriv15_I %*%
                                                                    solve(ImpCov) %*% SampCov))
  }



  #lapply(li,hess_fun(indexI=li$i,indexJ=li$j,Areg,A,ImpCov,SampCov,B,C,E,S))

 rett = matrix(unlist(snowfall::sfLapply(li,function(x) hess_fun(indexI=x["i"],indexJ=x["j"],Areg,A,ImpCov,SampCov,B,C,E,S))),length(par),length(par))
 #as.numeric(rett)
 rett
}
