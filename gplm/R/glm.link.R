"glm.link" <- function(eta, family="gaussian", link="identity", k=1){

  if (family=="bernoulli"){
    if (link=="logit"){
      mu<- 1/(1+exp(-eta))
      mu <- mu+ 1e-16*(mu==0) -1e-16*(mu==1)
    }
    if (link=="probit"){
      mu <- pnorm(eta)
      mu <- mu+ 1e-16*(mu==0) -1e-16*(mu==1)
    }
  }

  if (family=="gaussian"){
    if (link=="identity"){
      mu <- eta
    }
  }

  if (link=="log"){  ## gaussian+poisson
    mu <- exp(eta)
  }

  if (family=="gamma"){
    if (link=="reciprocal"){
      mu <- 1/eta
    }
  }

  if (family=="inverse.gaussian"){
    if (link=="quadratic.reciprocal"){
      mu <- 1/sqrt(eta)
    }
  }

  if (family=="negative.binomial"){
    if (link=="log"){
      e <- exp(eta)
      mu <- e/(k*(1-e))
    }
  }

  return(mu)
}

