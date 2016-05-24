theta.HB <-
function(result, subset){
     type <- result$type
     beta <- result$beta[,subset]
     sigv <- result$sigv[subset]
     Y <- result$Y
     X <- result$X
     m <- nrow(X)
     if (result$type == "FH"){
         Z <- result$Z
         b <- result$b
         phi <- outer(1 / b, sigv, "/")
         phi <- phi + 1 / Z
         phi <- 1 / (Z * phi)
         theta_HB <- phi * Y + (1 - phi) * (X %*% beta)
         theta_HB <- rowMeans(theta_HB)
     }
     else if (result$type == "YC"){
         sig2 <- result$sig2[,subset]
         b <- result$b
         phi <- outer(1 / b, sigv, "/")
         phi <- phi + 1 / sig2
         phi <- 1 / (sig2 * phi)
         theta_HB <- phi * Y + (1 - phi) * (X %*% beta)
         theta_HB <- rowMeans(theta_HB) 
     }
     else if (result$type == "SFH"){
         Z <- result$Z
         lam <- result$lam[subset]
         prox <- result$prox
         R <- array(0, c(m, m))
         R[prox] <- -1
         li1 <- prox[,1] 
         li2 <- prox[,2]
         R[cbind(li2, li1)] <- -1
         num <- rep(0, m)
         for (i in 1:m){
             num[i] = sum(li1 == i) + sum(li2 == i)
         }
         diag(R) <- num
         theta_HB = rowMeans(apply(rbind(beta, sigv, lam), MARGIN = 2, 
             FUN = theta.HBSFH, y = Y, X1 = X, R = R, vardir = Z))
     }
     else{
         sig2 <- result$sig2[,subset]
         lam <- result$lam[subset]
         prox <- result$prox
         R <- array(0, c(m, m))
         R[prox] <- -1
         li1 <- prox[,1] 
         li2 <- prox[,2]
         R[cbind(li2, li1)] <- -1
         num <- rep(0, m)
         for (i in 1:m){
             num[i] = sum(li1 == i) + sum(li2 == i)
         } 
         diag(R) <- num
         theta_HB = rowMeans(apply(rbind(beta, sigv, sig2, lam), MARGIN =  2, 
             FUN = theta.HBSYC, y = Y, X1 = X, R = R))
     }
     theta_HB
}
