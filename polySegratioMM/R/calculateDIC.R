`calculateDIC` <-
function(mcmc.mixture, model, priors, seg.ratios, chain=1, print.DIC=FALSE)
{

  ## calculate DIC_4 from Celeax et al
  
  if (class(mcmc.mixture) != "segratioMCMC")
    stop("'mcmc.mixture' must be of class 'segratioMCMC'")
  if (class(model) != "modelSegratioMM")
    stop("'model' must be of class 'modelSegratioMM'")
  if (class(priors) != "priorsSegratioMM")
    stop("'priors' must be of class 'priorsSegratioMM'")
  if (class(seg.ratios) != "segRatio") {
    stop("'seg.ratios' must be of class 'segRatio'")
  }

  mcmc <- mcmc.mixture$mcmc.list[[chain]]
  y <- gtools::logit(seg.ratios$seg.ratio)
  K <- model$n.components
  n <- length(y)
  m <- dim(mcmc)[1]
  
  p.ind <- grep("P\\[", colnames(mcmc))  # 1,2,3,4
  z.ind <- grep("T\\[", colnames(mcmc))  # 5 ...295
  mu.ind <- grep("mu\\[", colnames(mcmc))  #  296 297 298 299
  sigma.ind <- grep("sigma\\[", colnames(mcmc))  #  300 301 302 303
  
  equal.variances <- model$equal.variances   # modify if equal vars
  if (equal.variances)  {
    sigma.ind <- grep("sigma", colnames(mcmc))
  }

  ## need to think about this but it should be straight forward - loops
  ## to avoid massive memory
  
  ## I really need to learn to vectorise more of these calcs!
  
  first <- 0
  for (l in 1:m){ # l <- 1
    tij <- mcmc[l,p.ind]

    ynorm <- NULL
    if (equal.variances){
      for (j in 1:K){
        ynorm <- cbind(ynorm, pnorm(y, mcmc[l,mu.ind[j]],mcmc[l,sigma.ind]))
      }
    } else {
      for (j in 1:K){
        ynorm <- cbind(ynorm, pnorm(y, mcmc[l,mu.ind[j]],mcmc[l,sigma.ind[j]]))
      }
    }
###   for (j in 1:K){
###      ynorm <- cbind(ynorm, pnorm(y, mcmc[l,mu.ind[j]],mcmc[l,sigma.ind[j]]))
###    }
    
    tij <- mcmc[l,p.ind]*ynorm   
    tij <- tij/rowSums(tij)
    ## dim(mcmc[m,p.ind] %o% rep(1,n))
    first <- first + sum(tij* log( rep(1,n) %o% mcmc[l,p.ind] * ynorm),
                         na.rm=TRUE)
  }
  
  mu <- mcmc[,mu.ind]   # set up the original parameterisation
  for (i in length(mu.ind):2) {
    mu[,i] <- mu[,i] -  mu[,i-1] 
  }
  
  ## need priors for second sum on page 13 Celeaux et al (2006)
  mu.0 <- priors$params$logit.means
  var.0 <- 1/priors$params$logit.prec
  nu.0 <- 2*priors$params$A
  var.nu.0 <- 2* priors$params$B
  if (length(priors$params$dirichlet) == 0)
    alpha <- rep(1, length(p.ind))
  
  n.j <- rep(1,K)     # imagine that it should be 1 given specified priors
  
  mjl <- t(apply(mcmc[,z.ind], 1, table))
  p.jl <- (mjl+1)/(K+n)            # \bar{p}_j^(l)
  
  muhat.jl <- NA*mcmc[,p.ind]# calculate \hat{\mu_j^(l) and \bar{\mu_j^(l)
  s2hat.jl <- muhat.jl   # and also \hat{\sigma_j^2(l) and \bar{\sigma_j^2(l)
  for (i in 1:m){
    muhat.jl[i,] <- tapply(y, mcmc[i,z.ind], mean)
    for (j in 1:K){ # j <- 1
      s2hat.jl[i,j] <- sum((y[mcmc[i,z.ind]==j] - muhat.jl[i,j])^2)
    }
  }
  mu.jl <- ( rep(1,m) %o% (n.j*mu.0) + muhat.jl*mjl )/( rep(1,m) %o% n.j + mjl)
  if (equal.variances){
    s2.jl <- ( var.nu.0 + s2hat.jl  +
              (rep(1,m) %o% n.j * mjl)/(rep(1,m) %o% n.j + mjl)*
              (muhat.jl - rep(1,m)%o%mu.0)^2 /(nu.0 + mjl -2) )
  } else {
    s2.jl <- ((rep(1,m) %o% var.nu.0) + s2hat.jl  +
              (rep(1,m) %o% n.j * mjl)/(rep(1,m) %o% n.j + mjl)*
              (muhat.jl - rep(1,m)%o%mu.0)^2 /(rep(1,m)%o%nu.0 + mjl -2) )
  }
  
  second <- 0
  for (i in 1:m){
    second <- second + sum(log(p.jl[i,mcmc[i,z.ind]]) +   
                           pnorm(y, mu.jl[i,mcmc[i,z.ind]],
                                 s2.jl[i,mcmc[i,z.ind]], log.p = TRUE),
                           na.rm=TRUE)
  }
  
  DIC <- - (4/m) * first + (2/m)*second

  if (print.DIC)
    cat("DIC:",DIC,"\n")
  
  return(DIC)
}

