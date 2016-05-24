"glm.inverse.link" <- function(mu, family="gaussian", link="identity", k=1){

  if (family=="bernoulli"){
    if (link=="logit"){
      eta <- log(mu/(1-mu))
    }
    if (link=="probit"){
      eta <- qnorm(mu)
    }
  }

  if (family=="gaussian"){
    if (link=="identity"){
      eta <- mu
    }
  }

  if (link=="log"){  ## gaussian+poisson
    eta <- log(mu)
  }

  if (family=="gamma"){
    if (link=="reciprocal"){
      eta <- 1/mu
    }
  }

  if (family=="inverse.gaussian"){
    if (link=="quadratic.reciprocal"){
      eta <- 1/(mu*mu)
    }
  }

  if (family=="negative.binomial"){
    if (link=="log"){
      eta <- log(k*mu/(1+mu))
    }
  }

  return(eta)
}

