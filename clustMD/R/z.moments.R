z.moments <-
function(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, K, norms, nom.ind.Z){
    D <- J
    if(J > OrdIndx)
      D <- OrdIndx + sum(K[(OrdIndx+1):J] - 1)
    
    Ez.new <- array(NA, c(N, D, G))
    S.new <- matrix(0, N, G)
    S2.new <- array(NA, c(N, D, G))
    probs.new <- NA #dummy required for output
    if(J > OrdIndx)  probs.new <- array(NA, c(J-OrdIndx, max(K[(OrdIndx+1):J]), G))
    
    for(g in 1:G){
      #continuous
      if(CnsIndx > 0) Ez.new[, 1:CnsIndx, ] <- Y[, 1:CnsIndx]
      
      #ordinal
      if(OrdIndx > CnsIndx){
        for(i in 1:N){
          temp.e <- etruncnorm(a=zlimits[i, (CnsIndx+1):OrdIndx,1], b=zlimits[i, (CnsIndx+1):OrdIndx,2], mean=mu[(CnsIndx+1):OrdIndx, g], sd=sqrt(diag(Sigma[, , g])[(CnsIndx+1):OrdIndx]))
          Ez.new[i, (CnsIndx+1):OrdIndx, g] <- temp.e
          
          temp.v <- vtruncnorm(a=zlimits[i, (CnsIndx+1):OrdIndx,1], b=zlimits[i, (CnsIndx+1):OrdIndx,2], mean=mu[(CnsIndx+1):OrdIndx, g], sd=sqrt(diag(Sigma[, , g])[(CnsIndx+1):OrdIndx]))
          S.new[i, g] <- sum(temp.v) + t(temp.e)%*%temp.e
          S2.new[i, (CnsIndx+1):OrdIndx, g] <- temp.v + temp.e^2
        } # i
      } # if
      
      # Nominal
      if(J > OrdIndx){
        #       Nnorms <- nrow(norms)
        #       norms <- mvrnorm(Nnorms, mu=rep(0, max(K[(OrdIndx+1):J])-1), Sigma=diag(max(K[(OrdIndx+1):J])-1))
        for(j in (OrdIndx+1):J){
          Zrep <- norms[, 1:(K[j]-1)]%*%chol(Sigma[nom.ind.Z[[j-OrdIndx]], nom.ind.Z[[j-OrdIndx]], g]) + matrix(mu[nom.ind.Z[[j-OrdIndx]], g], dim(norms)[1], K[j]-1, byrow=TRUE)
          temp.z <- z.nom.diag(Zrep)
          probs.new[j-OrdIndx, 1:K[j], g] <- temp.z[[1]]
          for(k in 1:K[j]){
            Ez.new[Y[, j]==k, nom.ind.Z[[j-OrdIndx]], g] <- matrix(temp.z[[2]][, k], sum(Y[, j]==k), K[j]-1, byrow=TRUE)
            S.new[Y[, j]==k, g] <- S.new[Y[, j]==k, g] + sum(temp.z[[3]][, k])
            S2.new[Y[, j]==k, nom.ind.Z[[j-OrdIndx]], g] <- matrix(temp.z[[3]][, k], sum(Y[, j]==k), K[j]-1, byrow=TRUE)
          } # k
        } # j
      } # if
    } # g
    list(Ez.new, S.new, probs.new, S2.new)
  }
