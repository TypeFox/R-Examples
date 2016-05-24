bgmm.e.step <- function(X, model.params) {
  # densities for all components
  lfik <- matrix(0, nrow(X), model.params$k)
  for (i in 1:model.params$k) {
      if (model.params$d > 1) {
        ss = svd(model.params$cvar[i,,])
        rtas <- ss$d
        matc = t(ss$u[rtas > 10^-8, ]) %*% diag(rtas[rtas > 10^-8]^(-1/2)) %*% ss$v[rtas > 10^-8,]
        tx = apply(X, 1, get("-"), model.params$mu[i,,drop=F])
        lfik[,i] <-  -colSums((matc %*% tx)^2)/2 - sum(log(2*pi*rtas[rtas > 10^-8]))/2
      } else {
        lfik[,i] <-   dnorm(X, model.params$mu[i,,drop=F], sqrt(model.params$cvar[i,,]), log=T)
      }
  }
  # prior distribution 
  b.pi <- rbind(model.params$B, repeat.rows(model.params$pi, model.params$n-model.params$m))
  # normalisation step
  fb.ik = exp(lfik) * b.pi
  
  # numeric problems, trying to adjust, by keeping likelihood not so far from each other
  lcorrection <- 0
  if (max(fb.ik) == 0) {
    lcorrection <- max(lfik) 
    lfik = lfik - lcorrection
    # unbalanced case
    toCorr <- which(apply(lfik, 2, max) < -10)
    if (length(toCorr)>0) {
      for (correct in toCorr) 
        lfik[, correct] <- lfik[, correct] - max( lfik[, correct]) - 10
    }
    fb.ik = exp(lfik) * b.pi
  }
  # normalisation step
  
 list(tij =  t(apply(fb.ik, 1, normalize)), log.likelihood=sum(log(rowSums(fb.ik))) )
}

