bgmm.m.step <- function(X, model.params, model.structure, tij, eps=10^-6, priors.like.bgmm=TRUE) {
  new.model.params = model.params
 
  # new prior distribution
  if (priors.like.bgmm) {
    new.model.params$pi = colMeans(tij[-(1:model.params$m),])
  } else { # soft
    new.model.params$pi = colMeans(tij)
  }

  # new means 
  if  (model.structure$mean=="D") {
    # different mean for every component
    for (i in 1:model.params$k) 
         new.model.params$mu[i, ] = apply(X, 2, weighted.mean, tij[,i])
  } else {
    # same mean for every component
    new.model.params$mu = repeat.rows(colMeans(X), model.params$k)
  }

  # new variance matrix  
  for (i in 1:model.params$k) {
       tmp       = (X - repeat.rows(new.model.params$mu[i, ], model.params$n)) * sqrt(tij[,i])
       new.model.params$cvar[i, , ] = t(tmp) %*% tmp / sum(tij[,i])
       if (det(new.model.params$cvar[i, , ]) < eps)
          new.model.params$cvar[i, , ] = model.params$cvar[i, , ]
  }
  # are variances equal?
  if (model.structure$between=="E") {
      # averaging among clusters
     ncvar = matrix(0, model.params$d, model.params$d)
     for (i in 1:model.params$k) 
        ncvar = ncvar + new.model.params$cvar[i, , ] * new.model.params$pi[i]
     for (i in 1:model.params$k) 
        new.model.params$cvar[i, , ] = ncvar
  }
  if (model.structure$within=="E" && model.params$d>1) {
      # averaging among variables
     for (i in 1:model.params$k) {
        ndiag = sum(diag(new.model.params$cvar[i, , ]))
        sdiag = ndiag/model.params$d
        noutd = min(sdiag, (sum(new.model.params$cvar[i, , ])-ndiag)/(model.params$d*(model.params$d-1)))
        new.model.params$cvar[i, , ] = noutd
        diag(new.model.params$cvar[i, , ]) = sdiag
     }
  }
  # are covariance equal to 0?
  if (model.structure$cov=="0") {
   # covariance is equal to 0
   for (i in 1:model.params$k) 
      new.model.params$cvar[i, , ] = diag(diag(new.model.params$cvar[i, , ]), nrow=model.params$d)
  }
  
  new.model.params
}
