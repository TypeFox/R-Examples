
grad_ram_wMean = function(par,ImpCov,SampCov22,Areg,Sreg,
                           A,S,F,SampMean,lambda,type,m,mu,m.pars){


  #mats = extractMatrices(fit.growth)
  #A = mats$A
  #Areg = mats$A.est
  #S = mats$S
  #Sreg = mats$S.est
  #F = mats$F



  I = diag(nrow(Areg))
  #ImpCov = F %*% solve(I-Areg) %*% Sreg %*% t(solve(I-Areg)) %*% t(F)

  grad.out <- rep(0,length(par))

  B = solve(I - Areg)
  C = diag(nrow(ImpCov)) - (solve(ImpCov) %*% SampCov22)
  E = B %*% Sreg %*% t(B)
  dd=SampMean
  b = dd - mu

  #lik = log(((2*pi)^length(b)) * det(ImpCov)) + trace(solve(ImpCov) %*% SampCov22) +  t(b) %*% solve(ImpCov) %*% b
  #lik = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov))  - m


  A.iter <- max(A)


  #

  if(type=="none"){

    for(i in 1:length(grad.out)){

      A2 <- A == i;
      A2[A2==T] <- 1
      S2 <- S == i;
      S2[S2==T] <- 1
      m2 = m.pars == i
      m2[m2==T] <- 1

      deriv15 <- F %*% B %*% A2 %*% E %*% t(F) + F %*% B %*% S2 %*% t(B) %*% t(F)
      deriv20 <- F %*% B %*% A2 %*% B %*% m + F %*% B %*% m2

      grad.out[i]  <- trace(solve(ImpCov) %*% deriv15 %*% C) -
                      (t(b) %*% solve(ImpCov) %*% deriv15 + 2 * t(deriv20)) %*% solve(ImpCov) %*% b


    }

  }


  else if(type=="lasso"){
    for(i in 1:length(grad.out)){

      A2 <- A == i;
      A2[A2==T] <- 1
      S2 <- S == i;
      S2[S2==T] <- 1

      deriv15 <- F %*% B %*% A2 %*% E %*% t(F) + F %*% B %*% S2 %*% t(B) %*% t(F)
      # left out mean part
      grad.out[i]  <- trace(solve(ImpCov) %*% deriv15 %*% C) + if(i <= A.iter) lambda*sign(Areg[A==i]) else(0)# just penalize when A


    }

  }

  else if(type=="ridge"){
    for(i in 1:length(grad.out)){

      A2 <- A == i;
      A2[A2==T] <- 1
      S2 <- S == i;
      S2[S2==T] <- 1

      deriv15 <- F %*% B %*% A2 %*% E %*% t(F) + F %*% B %*% S2 %*% t(B) %*% t(F)
      # left out mean part
      grad.out[i]  <- trace(solve(ImpCov) %*% deriv15 %*% C) +
                      if(i <= A.iter) 2*lambda*Areg[A==i] else(0)


    }

  }




  grad.out
}
