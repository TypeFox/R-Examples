M.step <-
function(tau, N, Elz, J, OrdIndx, D, G, Y, CnsIndx, S2, model, a){
    # Mixing weights
    pi.vec.new <- apply(tau, 2, sum)/N
    
    # Mean vectors
    mu.new <- sweep(apply(Elz, c(2,3), sum), 2, apply(tau, 2, sum), "/")
    # Scale mean (ONLY REQUIRED FOR NOMINALS)
    if(J > OrdIndx) mu.new[(OrdIndx+1):D, ] <- sweep(as.matrix(mu.new[(OrdIndx+1):D, ]), 1, apply(sweep(as.matrix(mu.new[(OrdIndx+1):D, ]), 2, pi.vec.new, "*"), 1, sum), "-")
    #   mu.new <- mu.true
    
    # Sigma
    # Can probably vectorise this some way for speed
    if(model=="EII"){
      tausum <- apply(tau, 2, sum)
      Sigma.new <- array(0, c(D, D, G))
      lambda.new <-  matrix(1, G, D)
      
      # Lambda (Continuous & Ordinal)
      if(OrdIndx > 0){
        lambda.new_A <- -2*sum(diag(t(apply(array(Elz[, 1:OrdIndx, ], c(N, OrdIndx, G)), c(2,3), sum))%*%mu.new[1:OrdIndx, ])) + t(diag(t(mu.new[1:OrdIndx, ])%*%mu.new[1:OrdIndx, ]))%*%tausum
        if(CnsIndx > 0) lambda.new_A <- lambda.new_A + sum(diag(Y[, 1:CnsIndx]%*%t(Y[, 1:CnsIndx])))
        if(OrdIndx > CnsIndx) lambda.new_A <- lambda.new_A + sum(apply(array(S2[, (CnsIndx+1):OrdIndx, ], c(N, OrdIndx-CnsIndx, G)), c(1,3), sum)*tau)
        lambda.new_A <- as.numeric(lambda.new_A)/{N*OrdIndx}
        
        for(g in 1:G){
          Sigma.new[1:OrdIndx, 1:OrdIndx, g] <- lambda.new_A*diag(OrdIndx)
          lambda.new[g, 1:OrdIndx] <- lambda.new_A
        } 
      }
      
      if(J > OrdIndx){
        for(g in 1:G) Sigma.new[(OrdIndx+1):D, (OrdIndx+1):D, g] <- diag(D-OrdIndx)
      } 
      
      a.new <- matrix(1, G, D)
      
    }else if(model=="VII"){
      tausum <- apply(tau, 2, sum)
      Sigma.new <- array(0, c(D, D, G))
      lambda.new <-  matrix(1, G, D)
      
      if(OrdIndx > 0){
        lambda.new_A <- rep(NA, G)
        
        for(g in 1:G){
          lambda.new_A[g] <- -2*t(apply(Elz[,1:OrdIndx , g], 2, sum))%*%mu.new[1:OrdIndx, g] + t(mu.new[1:OrdIndx, g])%*%mu.new[1:OrdIndx, g]*tausum[g]
          if(CnsIndx > 0) lambda.new_A[g] <- lambda.new_A[g] + sum(diag(tau[, g]*Y[, 1:CnsIndx]%*%t(Y[, 1:CnsIndx])))
          if(OrdIndx > CnsIndx) lambda.new_A[g] <- lambda.new_A[g] + t(apply(S2[, (CnsIndx+1):OrdIndx, g], 1, sum))%*%tau[, g]
          lambda.new_A[g] <- lambda.new_A[g]/{tausum[g]*OrdIndx}
          lambda.new[g, 1:OrdIndx] <- lambda.new_A[g]
          Sigma.new[1:OrdIndx, 1:OrdIndx, g] <- lambda.new_A[g]*diag(OrdIndx)
        } # g 
      }
      
      if(J > OrdIndx){
        lambda.new_B <- rep(NA, G)
        for(g in 1:G){
          lambda.new_B[g] <- -2*t(apply(Elz[, (OrdIndx+1):D, g], 2, sum))%*%mu.new[(OrdIndx+1):D, g] + t(mu.new[(OrdIndx+1):D, g])%*%mu.new[(OrdIndx+1):D, g]*tausum[g] + t(apply(S2[, (OrdIndx+1):D, g], 1, sum))%*%tau[, g]
          lambda.new_B[g] <- lambda.new_B[g]/{tausum[g]*(D-OrdIndx)}
        } # g
        lambda.new_B <- lambda.new_B/sum(lambda.new_B)
        
        for(g in 1:G){
          lambda.new[g, (OrdIndx+1):D] <- lambda.new_B[g]
          Sigma.new[(OrdIndx+1):D, (OrdIndx+1):D, g] <- lambda.new_B[g]*diag(D-OrdIndx)
        } # g
      }
      
      a.new <- matrix(1, G, D)
      
    }else if(model=="EEI"){
      
      tausum <- apply(tau, 2, sum)
      Sigma.new <- array(0, c(D, D, G))
      lambda.new <-  matrix(1, G, D)
      
      # lambda (continuous & ordinal)
      if(OrdIndx > 0){
        lambda.new_A <- -2*sum(diag(t(apply(array(Elz[, 1:OrdIndx, ], c(N, OrdIndx, G)), c(2,3), sum))%*%diag(1/a[1,1:OrdIndx])%*%mu.new[1:OrdIndx, ])) + sum(diag(t(mu.new[1:OrdIndx, ])%*%diag(1/a[1, 1:OrdIndx])%*%mu.new[1:OrdIndx, ])*tausum)
        if(CnsIndx > 0) lambda.new_A <- lambda.new_A + sum(sweep(matrix(Y[, 1:CnsIndx]^2, ncol=CnsIndx),2,a[1,1:CnsIndx],"/"))
        if(OrdIndx > CnsIndx) lambda.new_A <- lambda.new_A + sum(apply(sweep(array(S2[, (CnsIndx+1):OrdIndx, ], c(N, OrdIndx-CnsIndx, G)), c(1,3), tau, "*"), c(1,2), sum)%*%(1/a[1, (CnsIndx+1):OrdIndx]))
        lambda.new_A <- lambda.new_A/{N*OrdIndx}
        lambda.new[, 1:OrdIndx] <- lambda.new_A
      }
      
      # A
      a.new <- -2*diag(apply(Elz, c(2, 3), sum)%*%t(mu.new)) + mu.new^2%*%tausum
      if(CnsIndx > 0) a.new[1:CnsIndx] <- a.new[1:CnsIndx] + apply(matrix(Y[, 1:CnsIndx]^2, N, ncol=CnsIndx), 2, sum)
      if(J > CnsIndx) a.new[(CnsIndx+1):D] <- a.new[(CnsIndx+1):D] + apply(apply(sweep(array(S2[, (CnsIndx+1):D, ], c(N, D-CnsIndx, G)), c(1,3), tau, "*"), c(1,2), sum), 2, sum)
      a.new <- a.new/(N*lambda.new[1, ])
      # det(A) = 1 constraint
      if(OrdIndx > 0)
        a.new[1:OrdIndx] <- a.new[1:OrdIndx]/prod(a.new[1:OrdIndx])^(1/OrdIndx)
      
      if(J > OrdIndx) a.new[(OrdIndx+1):D] <- 1
      
      a.new <- matrix(a.new, G, D, byrow=TRUE)
      Sigma.new <- array(diag(lambda.new[1, ]*a.new[1, ]), c(D, D, G))
      
    }else if(model=="VEI"){
      tausum <- apply(tau, 2, sum)
      Sigma.new <- array(0, c(D, D, G))
      lambda.new <-  matrix(1, G, D)
      
      # lambda (continuous & ordinal)
      if(OrdIndx > 0){
        lambda.new_A <- -2*diag(t(apply(array(Elz[, 1:OrdIndx, ], c(N, OrdIndx, G)), c(2, 3), sum))%*%diag(1/a[1, 1:OrdIndx])%*%mu.new[1:OrdIndx, ]) + (t(1/a[1, 1:OrdIndx])%*%{mu.new[1:OrdIndx, ]^2})*tausum
        if(CnsIndx > 0) lambda.new_A <- lambda.new_A + t((Y[, 1:CnsIndx]^2)%*%(1/a[1, 1:CnsIndx]))%*%tau
        if(OrdIndx > CnsIndx) lambda.new_A <- lambda.new_A + (1/a[1, (CnsIndx+1):OrdIndx])%*%apply(sweep(array(S2[, (CnsIndx+1):OrdIndx, ], c(N, OrdIndx-CnsIndx, G)), c(1,3), tau, "*"), c(2, 3), sum)
        lambda.new[, 1:OrdIndx] <- lambda.new_A/{OrdIndx*tausum}
      }
      
      # lambda (nominal)
      if(J > OrdIndx){
        lambda.new_B <- -2*diag(t(apply(array(Elz[, (OrdIndx+1):D, ], c(N, D-OrdIndx, G)), c(2, 3), sum))%*%diag(1/a[1, (OrdIndx+1):D])%*%mu.new[(OrdIndx+1):D, ]) + (t(1/a[1, (OrdIndx+1):D])%*%{mu.new[(OrdIndx+1):D, ]^2})*tausum + (1/a[1, (OrdIndx+1):D])%*%apply(sweep(array(S2[, (OrdIndx+1):D, ], c(N, D-OrdIndx, G)), c(1,3), tau, "*"), c(2, 3), sum)
        lambda.new_B <- lambda.new_B/{(D-OrdIndx)*tausum}
        lambda.new[, (OrdIndx+1):D] <- lambda.new_B/sum(lambda.new_B)
      }
      
      # A
      a.new <- -2*diag(apply(Elz, c(2, 3), sum)%*%t(mu.new/t(lambda.new))) + t({mu.new^2/t(lambda.new)}%*%tausum)
      if(CnsIndx > 0) a.new[1:CnsIndx] <- a.new[1:CnsIndx] + (1/lambda.new[, 1])%*%t(tau)%*%(Y[, 1:CnsIndx]^2)
      if(J > CnsIndx) a.new[(CnsIndx+1):D] <- a.new[(CnsIndx+1):D] + diag(apply(sweep(array(S2[, (CnsIndx+1):D, ], c(N, D-CnsIndx, G)), c(1,3), tau, "*"), c(2, 3), sum)%*%(1/lambda.new[, (CnsIndx+1):D]))
      a.new <- a.new/N
      # det(A) = 1 constraint
      if(OrdIndx > 0)
        a.new[1:OrdIndx] <- a.new[1:OrdIndx]/prod(a.new[1:OrdIndx])^(1/OrdIndx)
      
      # identifiability constraint
      if(J > OrdIndx)  a.new[(OrdIndx+1):D] <- 1
      
      a.new <- matrix(a.new, G, D, byrow=TRUE)
      Sigma.new <- array(NA, c(D, D, G))
      for(g in 1:G)
        Sigma.new[, , g] <- lambda.new[g]*diag(a.new[1, ])
      
    }else if(model=="EVI"){
      tausum <- apply(tau, 2, sum)
      Sigma.new <- array(0, c(D, D, G))
      lambda.new <-  matrix(1, G, D)
      
      # lambda (continuous & ordinal)
      if(OrdIndx > 0){
        lambda.new_A <- -2*sum(apply(array(Elz[, 1:OrdIndx, ], c(N, OrdIndx, G)), c(2, 3), sum)*matrix(mu.new[1:OrdIndx, ]/t(a[, 1:OrdIndx]), OrdIndx, G)) + diag(t(mu.new[1:OrdIndx, ]^2)%*%(1/t(matrix(a[, 1:OrdIndx], G, OrdIndx))))%*%tausum
        if(CnsIndx > 0) lambda.new_A <- lambda.new_A + sum(diag(t(tau)%*%(Y[, 1:CnsIndx]^2)%*%(t(matrix(1/a[, 1:CnsIndx], G, CnsIndx)))))
        if(OrdIndx > CnsIndx) lambda.new_A <- lambda.new_A + sum(diag(t(apply(sweep(array(S2[, (CnsIndx+1):OrdIndx, ], c(N, OrdIndx-CnsIndx, G)), c(1, 3), tau, "*"), c(2, 3), sum))%*%t(matrix(1/a[, (CnsIndx+1):OrdIndx], G, OrdIndx-CnsIndx))))
        lambda.new_A <- lambda.new_A/{OrdIndx*N}
        lambda.new[, 1:OrdIndx] <- lambda.new_A
      }
      
      # A
      a.new <- t(-2*apply(Elz, c(2,3), sum)*mu.new + sweep(mu.new^2, 2, tausum, "*"))
      if(CnsIndx > 0) a.new[, 1:CnsIndx] <- a.new[, 1:CnsIndx] + t(tau)%*%(Y[, 1:CnsIndx]^2)
      if(J > CnsIndx) a.new[, (CnsIndx+1):D] <- a.new[, (CnsIndx+1):D] + t(apply(sweep(array(S2[, (CnsIndx+1):D, ], c(N, D-CnsIndx, G)), c(1, 3), tau, "*"), c(2, 3), sum))
      a.new <- a.new/sweep(lambda.new, 1, tausum, "*")
      # det(A) = 1 constraint
      if(OrdIndx > 0)
        a.new[, 1:OrdIndx] <- sweep(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, apply(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, prod)^(1/OrdIndx), "/")
      
      #scale
      if(J > OrdIndx) a.new[, (OrdIndx+1):D] <- sweep(matrix(a.new[, (OrdIndx+1):D], G, D-OrdIndx), 2, apply(matrix(a.new[, (OrdIndx+1):D], G, D-OrdIndx), 2, sum), "/")
      
      Sigma.new <- array(NA, c(D, D, G))
      for(g in 1:G)
        Sigma.new[, , g] <- lambda.new[g, ]*diag(a.new[g, ])
      
    }else if(model=="VVI"){
      tausum <- apply(tau, 2, sum)
      Sigma.new <- array(0, c(D, D, G))
      lambda.new <-  matrix(1, G, D)
      
      # lambda (continuous & ordinal)
      if(OrdIndx > 0){
        lambda.new_A <- -2*diag(t(apply(array(Elz[, 1:OrdIndx, ], c(N, OrdIndx, G)), c(2,3), sum))%*%matrix(mu.new[1:OrdIndx, ]/t(a[, 1:OrdIndx]), OrdIndx, G)) + diag(t(mu.new[1:OrdIndx, ]^2)%*%(1/t(matrix(a[, 1:OrdIndx], G, OrdIndx))))*tausum
        if(CnsIndx > 0) lambda.new_A <- lambda.new_A + diag(t(tau)%*%(Y[, 1:CnsIndx]^2)%*%(t(matrix(1/a[, 1:CnsIndx], G, CnsIndx))))
        if(OrdIndx > CnsIndx) lambda.new_A <- lambda.new_A + diag(t(apply(sweep(array(S2[, (CnsIndx+1):OrdIndx, ], c(N, OrdIndx-CnsIndx, G)), c(1, 3), tau, "*"), c(2, 3), sum))%*%t(matrix(1/a[, (CnsIndx+1):OrdIndx], G, OrdIndx-CnsIndx)))
        lambda.new_A <- lambda.new_A/{OrdIndx*tausum}
        lambda.new[, 1:OrdIndx] <- lambda.new_A      
      }
      
      # lambda (nominal)
      if(J > OrdIndx){
        lambda.new_B <- -2*diag(t(apply(array(Elz[, (OrdIndx+1):D, ], c(N, D-OrdIndx, G)), c(2,3), sum))%*%matrix(mu.new[(OrdIndx+1):D, ]/t(a[, (OrdIndx+1):D]), D-OrdIndx, G)) + diag(t(mu.new[(OrdIndx+1):D, ]^2)%*%(1/t(matrix(a[, (OrdIndx+1):D], G, D-OrdIndx))))*tausum + diag(t(apply(sweep(array(S2[, (OrdIndx+1):D, ], c(N, D-OrdIndx, G)), c(1, 3), tau, "*"), c(2, 3), sum))%*%t(matrix(1/a[, (OrdIndx+1):D], G, D-OrdIndx)))
        lambda.new_B <- lambda.new_B/{(D-OrdIndx)*tausum}
        lambda.new[, (OrdIndx+1):D] <- lambda.new_B/sum(lambda.new_B)
      }
      
      # A
      a.new <- t(-2*apply(Elz, c(2, 3), sum)*mu.new + sweep(mu.new^2, 2, tausum, "*"))
      if(CnsIndx > 0) a.new[, 1:CnsIndx] <- a.new[, 1:CnsIndx] + t(tau)%*%(Y[, 1:CnsIndx]^2)
      if(J > CnsIndx) a.new[, (CnsIndx+1):D] <- a.new[, (CnsIndx+1):D] + t(apply(sweep(array(S2[, (CnsIndx+1):D, ], c(N, D-CnsIndx, G)), c(1, 3), tau, "*"), c(2, 3), sum))
      a.new <- a.new/sweep(lambda.new, 1, tausum, "*")
      # det(A) = 1 constraint
      if(OrdIndx > 0)
        a.new[, 1:OrdIndx] <- sweep(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, apply(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, prod)^(1/OrdIndx), "/")
      
      #scale
      if(J > OrdIndx) a.new[, (OrdIndx+1):D] <- sweep(matrix(a.new[, (OrdIndx+1):D], G, D-OrdIndx), 2, apply(matrix(a.new[, (OrdIndx+1):D], G, D-OrdIndx), 2, sum), "/")
      
      for(g in 1:G)
        Sigma.new[, , g] <- lambda.new[g, ]*diag(a.new[g, ])
      
    }else{
      
      print("Unknown model chosen! Choose from; EII, VII, EEI, VEI, EVI or VVI")
      
    }  
    
    list(pi.vec.new, mu.new, lambda.new, a.new, Sigma.new)
  }
