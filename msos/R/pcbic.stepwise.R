pcbic.stepwise <-
function(eigenvals,n) {
     k <- length(eigenvals)
     p0 <- rep(1,k)
     b <- rep(0,k)
     pb <- vector('list',k)
     pb[[1]] <- p0
     b[1] <- pcbic(eigenvals,n,p0)$BIC
     for(i in 2:k) {
          psb <- pcbic.subpatterns(eigenvals,n,pb[[i-1]])
          b[i] <- min(psb$bic)
          pb[[i]] <- psb$pattern[,psb$bic==b[i]]
     }
     ib <- (1:k)[b==min(b)]
     list(Patterns = pb,BICs = b,
           BestBIC = b[ib],BestPattern = pb[[ib]],
               LambdaHat = pcbic(eigenvals,n,pb[[ib]])$lambdaHat)
}
