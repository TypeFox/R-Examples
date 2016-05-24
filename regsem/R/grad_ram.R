
grad_ram = function(par,ImpCov,SampCov,Areg,Sreg,A,S,
                     F,lambda,type,pars_pen,diff_par){

  grad.out <- rep(0,length(par))

  B = solve(diag(nrow(A)) - Areg)
  C = diag(nrow(ImpCov)) - solve(ImpCov) %*% SampCov
  E = B %*% Sreg %*% t(B)



  # The S matrix gradients are exactly twice that of other methods

  if(type=="none"){

    for(i in 1:length(grad.out)){

      A2 <- A == i;
      A2[A2==T] <- 1
      S2 <- S == i;
      S2[S2==T] <- 1

      deriv15 <- F %*% B %*% A2 %*% E %*% t(F) + F %*% B %*% S2 %*% t(B) %*% t(F)
      # left out mean part
      grad.out[i]  <- trace(solve(ImpCov) %*% deriv15 %*% C)


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
      grad.out[i]  <- trace(solve(ImpCov) %*% deriv15 %*% C) + if(any(i==pars_pen)) lambda*sign(Areg[A==i]) else(0)# just penalize when A


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
                      if(any(i==pars_pen)) 2*lambda*Areg[A==i] else(0)


    }

  }

  else if(type=="diff_lasso"){
    count=0
    for(i in 1:length(grad.out)){

      A2 <- A == i;
      A2[A2==T] <- 1
      S2 <- S == i;
      S2[S2==T] <- 1

      deriv15 <- F %*% B %*% A2 %*% E %*% t(F) + F %*% B %*% S2 %*% t(B) %*% t(F)
      # left out mean part
      grad.out[i]  <- trace(solve(ImpCov) %*% deriv15 %*% C) +
        if(any(i==pars_pen)){
          count=count+1
          lambda*sign(Areg[A==i]-diff_par[count])
        }else(0)


    }

  }



  grad.out[(max(A)+1):max(S)] = grad.out[(max(A)+1):max(S)] *0.5
  #grad.out[min(S[S!=0],0):max(S)] = grad.out[min(S[S!=0],0):max(S)] *0.5
  grad.out
}
