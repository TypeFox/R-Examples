ObsLogLikelihood <-
function(N, CnsIndx, G, Y, mu, Sigma, pi.vec, J, OrdIndx, K, perc.cut, Nnorms, zlimits, nom.ind.Z){
    #Continuous
    LikeCns <- rep(1, N)
    if(CnsIndx > 0){
      densCns <- matrix(NA, N, G)
      for(g in 1:G){
        densCns[, g] <- dmvnorm(Y[, 1:CnsIndx], mean=mu[1:CnsIndx, g], sigma=matrix(Sigma[1:CnsIndx, 1:CnsIndx, g], CnsIndx, CnsIndx) )
      }
      densCns <- sweep(densCns, 2, pi.vec, "*")
      LikeCns <- apply(densCns, 1, sum)
    }
    
    # Categorical
    LikeCat <- rep(1, N)
    if(J > CnsIndx){
      # Ordinal
      densOrd <- matrix(1, N, G)
      if(OrdIndx > CnsIndx){
        OrdProbs <- array(NA, c(OrdIndx-CnsIndx, max(K[(CnsIndx+1):OrdIndx]),G))
        for(j in (CnsIndx+1):OrdIndx){
          for(g in 1:G){
            CumulProbs <- pnorm(perc.cut[[j]], mean=mu[j, g], sd=sqrt(Sigma[j, j, g]))
            for(k in 1:K[j])
              OrdProbs[j-CnsIndx, k, g] <- CumulProbs[k+1] - CumulProbs[k]
            
            densOrd[, g] <- densOrd[, g]*OrdProbs[j-CnsIndx, Y[, j], g]
          } # g
        } # j
      }
      
      # Nominal
      densNom <- matrix(1, N, G)
      if(J > OrdIndx){
        norms <- mvrnorm(Nnorms, mu=rep(0, max(K[(OrdIndx+1):J])-1), Sigma=diag(max(K[(OrdIndx+1):J])-1))
        # Z moments
        temp.z <- z.moments(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, K, norms, nom.ind.Z)
        probs.nom_obsLike <- temp.z[[3]]
        for(g in 1:G){
          for(j in (OrdIndx+1):J)
            densNom[, g] <- densNom[, g]*probs.nom_obsLike[j-OrdIndx, Y[, j], g]
        } # g
      }
      
      densCat <- densOrd*densNom
      densCat <- sweep(densCat, 2, pi.vec, "*")
      LikeCat <- apply(densCat, 1, sum)
    }
    logLike <- sum(log(LikeCns) + log(LikeCat))
    logLike
  }
