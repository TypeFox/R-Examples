E.step <-
function(N, G, pi.vec, Y, OrdIndx, CnsIndx, D, perc.cut, mu, Sigma, Ez, J, probs.nom, K){
    # tau
    tau.mat <- matrix(NA, N, G)
    
    for(g in 1:G){
      # continuous
      if(CnsIndx > 0){
        tau.cns <- dmvnorm(matrix(Y[, 1:CnsIndx], nrow=N), mu[1:CnsIndx, g], matrix(Sigma[1:CnsIndx, 1:CnsIndx, g], nrow=CnsIndx))
      }else{
        tau.cns <- rep(1, N)
      } 
      
      #ordinal
      if(OrdIndx > CnsIndx){
        tau.ord <- matrix(NA, N, OrdIndx-CnsIndx)
        for(j in (CnsIndx+1):OrdIndx){
          resp.probs.g <- pnorm(perc.cut[[j]][2:length(perc.cut[[j]])], mean=mu[j, g], sd=sqrt(Sigma[j, j, g])) - pnorm(perc.cut[[j]][1:(length(perc.cut[[j]])-1)], mean=mu[j, g], sd=sqrt(Sigma[j, j, g]))
          tau.ord[, j-CnsIndx] <- resp.probs.g[Y[,j]]
        }
        tau.ord <- apply(tau.ord, 1, prod)
      }else{
        tau.ord <- rep(1, N)
      }
      
      # nominal
      if(J > OrdIndx){
        tau.nom <- matrix(NA, N, J-OrdIndx)
        for(j in (OrdIndx+1):J){
          resp.probs.g <- probs.nom[j-OrdIndx, 1:K[j], g]
          tau.nom[, j-OrdIndx] <- resp.probs.g[Y[,j]]
        }
        tau.nom <- apply(tau.nom, 1, prod)
      }else{
        tau.nom <- rep(1, N)
      }
      
      tau.mat[, g] <- pi.vec[g]*tau.cns*tau.ord*tau.nom
      # should change function to work on log scale
    } # g
    tau.new <- sweep(tau.mat, 1, apply(tau.mat, 1, sum), "/")
    
    # E(lz)
    Elz.new <- sweep(Ez, c(1,3) , tau.new, "*")
    
    list(tau.new, Elz.new)
  }
