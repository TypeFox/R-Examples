epi.cpresids <- function(obs, fit, covpattern){
   # Covariate pattern identifiers:
   cpid <- covpattern$cov.pattern$id
   
   # Number of observations that comprise each covariate pattern:
   n <- covpattern$cov.pattern$n
   
   # Number of outcome-positive observations observed for each covariate pattern:
   nY <- obs

   # Number of outcome-positive observations predicted for each covariate pattern:
   np <- fit * n
   
   # Predicted probability of outcome for each covariate pattern:
   pi. <- fit

   den <- rep(1, times = nrow(covpattern$cov.pattern))

   # Turn factors into dummy variables:
   X <- den
   for(i in 3:dim(covpattern$cov.pattern)[2]){
   ifelse(is.factor(covpattern$cov.pattern[,i]), 
      # The function model.matrix returns the dummy variables for each factor. Remove the first column to return treatment contrasts.
      # That is, if you have a factor comprised of three levels, we return two columns to represent the treatment contrasts (i.e. 00, 01, and 10).
      X <- cbind(X, model.matrix(~covpattern$cov.pattern[,i] - 1)[,-1]),
      X <- cbind(X, covpattern$cov.pattern[,i]))
   }   
   colnames(X) <- 1:dim(X)[2]
   # X <- as.matrix(cbind(den, covpattern$cov.pattern[3:dim(covpattern$cov.pattern)[2]]))
   V <- diag(np * (1 - pi.))
   xvx <- solve(t(X) %*% V %*% X)
   sV <- sqrt(V)
   H <- sV %*% X %*% xvx %*% t(X) * sV
   leverage <- diag(H)

   # Raw residuals:
   raw <- (nY - np)
   
   # Standardised raw residuals:
   sraw <- raw /sd(np)
   
   # Pearson residuals:
   pearson <- (nY - np)/sqrt(np * (1 - pi.))

   # Standardised Pearson residuals:
   spearson <- pearson / sqrt(1 - leverage)

   # Deviance residuals:
   sign <- ifelse(nY - np > 0, 1, -1)
   dev <- sign * sqrt(2 * ((nY * log(nY/np)) + ((n - nY) * log((n - nY)/(n * (1 - pi.))))))
   dev[nY == 0] <- -sqrt(2 * n[nY == 0] * abs(log(1 - pi.[nY == 0])))
   dev[nY == n] <- sqrt(2 * n[nY == n] * abs(log(pi.[nY == n])))
   
   # Delta beta:
   deltabeta <- (pearson^2 * leverage) / (1 - leverage)

   # Standardised delta beta:
   sdeltabeta <- (spearson^2 * leverage) / (1 - leverage)
   
   # Delta chi-square (used to detect ill-fitting covariate patterns):
   deltachi <- pearson^2 / (1 - leverage)

   rval <- data.frame(cpid = cpid, n = n, obs = nY, pred = np, 
      raw = raw, sraw = sraw, pearson = pearson, spearson = spearson, 
      deviance = dev, leverage = leverage, deltabeta = deltabeta, sdeltabeta = sdeltabeta, deltachi = deltachi)
   return(rval)
}


