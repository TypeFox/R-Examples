soft.e.step <- function(X, model.params) {
  # densities for all components
#  fik  <- t(apply(X, 1, get.all.densities, model.params, model.params$k))
  lfik <- matrix(0, nrow(X), model.params$k)
  for (i in 1:model.params$k) {
      if (model.params$d > 1) {
        ss = svd(model.params$cvar[i,,])
        rtas <- ss$d
        matc = t(ss$u[rtas > 10^-8, ]) %*% diag(rtas[rtas > 10^-8]^(-1/2)) %*% ss$v[rtas > 10^-8,]
#        matc = t(ss$u) %*% diag(ss$d^(-1/2)) %*% ss$v
        tx = apply(X, 1, get("-"), model.params$mu[i,,drop=F])
#        lfik[,i] <-  -colSums((matc %*% tx)^2)/2 - sum(log(2*pi*ss$d))/2
        lfik[,i] <-  -colSums((matc %*% tx)^2)/2 - sum(log(2*pi*rtas[rtas > 10^-8]))/2
      } else {
#        tx = apply(X, 1, get("-"), model.params$mu[i,,drop=F])
#        lfik[,i] <-  -(tx^2)/(2*model.params$cvar[i,,]) - log(2*pi*model.params$cvar[i,,])/2
        lfik[,i] <-   dnorm(X, model.params$mu[i,,drop=F], sqrt(model.params$cvar[i,,]), log=T)
      }
  }  
        
  if (is.null(model.params$P)) { #unsupervised
    p.ik <- matrix(1,model.params$n,model.params$k)
  } else { #soft
    p.ik <- rbind(model.params$P, matrix(1/model.params$k,model.params$n-model.params$m,model.params$k))
  }
#  repeat.rows(rep(-log(model.params$k), model.params$k), model.params$n-model.params$m))

  fik = exp(lfik)
  # numeric problems, trying to adjust, by keeping likelihood not so far from each other
  lcorrection <- 0
  if (max(fik) == 0) {
  	lcorrection <- max(lfik) 
  	lfik = lfik - lcorrection
  	# unbalanced case
  	toCorr <- which(apply(lfik, 2, max) < -10)
  	if (length(toCorr)>0) {
  		for (correct in toCorr) 
  		  lfik[, correct] <- lfik[, correct] - max( lfik[, correct]) - 10
  	}
    fik = exp(lfik)
  }
  tij = t(t(p.ik) * model.params$pi) * fik
  margintij <- pmax(rowSums(tij), 10*.Machine$double.eps) # machine epsilon
  n.ltij = log(margintij)
  log.likelihood = sum(n.ltij) #  + lcorrection
  tij = tij/exp(n.ltij)
  # normalisation step
  
  list(tij =  tij, log.likelihood=log.likelihood )
}
