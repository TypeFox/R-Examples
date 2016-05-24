
hess_ram = function(par,ImpCov,SampCov,Areg,Sreg,A,S,F){

#?sfClusterApplyLB
#sfClusterApplyLB(index, simulationIteration, n, rho, g, h, pv, het=F)

hess.out <- matrix(0,length(par),length(par))
h <- 0.00001


B = solve(diag(nrow(A)) - Areg)
C = diag(nrow(ImpCov)) - solve(ImpCov) %*% SampCov
E = B %*% Sreg %*% t(B)

#if(type=="none"){

# not symmetric

    for(i in 1:nrow(hess.out)){
          for(j in 1:ncol(hess.out)){

          Ai <- (A == i)*1
          Aj <- (A == j)*1
          Aij <- matrix(0,nrow(A),ncol(A))
          Aij[Ai==T & Aj ==T] <- 1

          Si <- (S == i)*1;
          Sj <- (S == j)*1;
          Sij <- matrix(0,nrow(S),ncol(S))
          Sij[Si==T & Sj==T] <- 1

  deriv15_I <- F %*% B %*% Ai %*% E %*% t(F) + F %*% B %*% Si %*% t(B) %*% t(F)
  deriv15_J <- F %*% B %*% Aj %*% E %*% t(F) + F %*% B %*% Sj %*% t(B) %*% t(F)

  # this is cause of asymmetry
  deriv15_IJ <- F %*% B %*% Ai %*% B %*% Aj %*% E %*% t(F) +
                F %*% B %*% Aj %*% B %*% Ai %*% E %*% t(F) +
                F %*% B %*% Ai %*% B %*% Sj %*% t(B) %*% t(F) +
                F %*% B %*% Ai %*% E %*% t(Aj) %*% t(B) %*% t(F) +
                F %*% B %*% Aj %*% B %*% Si %*% t(B) %*% t(F)

  # left out mean part
  hess.out[i,j]  <- trace(solve(ImpCov) %*% deriv15_IJ %*% C - solve(ImpCov) %*%
                          deriv15_I %*% solve(ImpCov) %*% deriv15_J %*% C +
                          solve(ImpCov) %*% deriv15_J %*% solve(ImpCov) %*% deriv15_I %*%
                          solve(ImpCov) %*% SampCov)

          }
      }


hess.out



}
