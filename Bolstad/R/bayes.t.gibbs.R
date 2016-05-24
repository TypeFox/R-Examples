#' @importFrom utils tail
bayes.t.gibbs = function(x, y, nIter = 10000, nBurn = 1000, sigmaPrior = c("chisq", "gamma")){
  
  sigmaPrior = match.arg(sigmaPrior)
  nx = length(x)
  ny = length(y)
  xbar = mean(x)
  ybar = mean(y)
  Sx = nx * xbar
  Sy = ny * ybar
  SSx = sum((x - xbar)^2)
  SSy = sum((y - ybar)^2)
  
  if(sigmaPrior == "chisq"){
    ## prior mean
    m0x = median(x)
    m0y = median(y)
    
    ## prior sd
    s0x = sd(x)
    s0y = sd(y)
    
    ##
    S0x = s0x^2 * qchisq(0.5, 1)
    S0y = s0y^2 * qchisq(0.5, 1)
    
    ## 
    kappa = 1
    kappa1x = kappa + nx
    kappa1y = kappa + ny
    
    N = nIter + nBurn
    
    sigma.sq.x = sigma.sq.y = mu.x = mu.y = rep(0, N)
    
    ## draw the initial values
    sigma.sq.x[1] = S0x / rchisq(1, kappa1x)
    sigma.sq.y[1] = S0y / rchisq(1, kappa1y)
    mu.x[1] = rnorm(1, m0x, s0x)
    mu.y[1] = rnorm(1, m0y, s0y)
    
    
    for(i in 2:N){
      S1x = S0x + sum((x - mu.x[i - 1])^2)
      S1y = S0y + sum((y - mu.y[i - 1])^2)
      
      sigma.sq.x[i] = S1x / rchisq(1, kappa1x)
      sigma.sq.y[i] = S1y / rchisq(1, kappa1y)
      
      s1x = 1 / sqrt(1 / s0x^2 + nx / sigma.sq.x[i]) 
      s1y = 1 / sqrt(1 / s0y^2 + ny / sigma.sq.y[i])
      
      m1x = m0x * s1x^2  / s0x^2  + Sx * s1x^2 / sigma.sq.x[i] 
      m1y= m0y * s1y^2  / s0y^2  + Sy * s1y^2 / sigma.sq.y[i]
      
      mu.x[i] = rnorm(1, m1x, s1x) 
      mu.y[i] = rnorm(1, m1y, s1y) 
    }
    
    res = data.frame(mu.x = tail(mu.x, nIter),
                     mu.y = tail(mu.y, nIter),
                     mu.diff = tail(mu.x - mu.y, nIter),
                     sigma.sq.x = tail(sigma.sq.x, nIter),
                     sigma.sq.y = tail(sigma.sq.y, nIter),
                     tstat = tail((mu.x - mu.y) / sqrt(sigma.sq.x / nx + sigma.sq.y / ny), nIter))
    
  }else{
    ## prior means and sds
    ## prior mean
    m0x = median(x)
    m0y = median(y)
    
    ## prior sd
    s0x = sd(x)
    s0y = sd(y)
    
    alpha0x = beta0x = alpha0y = beta0y = 0.001
    
    alpha1x = nx / 2 + alpha0x
    beta1x = beta0x + SSx / 2
    
    alpha1y = ny / 2 + alpha0y
    beta1y = beta0y + SSy / 2
    
    N = nIter + nBurn
    
    sigma.sq.x = sigma.sq.y = mu.x = mu.y = rep(0, N)
    
    ## draw the initial values
    sigma.sq.x[1] = 1/rgamma(1, alpha1x, beta1x)
    sigma.sq.y[1] = 1/rgamma(1, alpha1y, beta1y)
    
    mu.x[1] = rnorm(1, m0x, s0x)
    mu.y[1] = rnorm(1, m0y, s0y)

    for(i in 2:N){
      beta1x = beta0x + sum((x - mu.x[i - 1])^2) * 0.5
      sigma.sq.x[i] = 1 / rgamma(1, alpha1x, beta1x)
      
      beta1y = beta0y + sum((y - mu.y[i - 1])^2) * 0.5
      sigma.sq.y[i] = 1 / rgamma(1, alpha1y, beta1y)
      
      s1x = 1 / sqrt(1 / s0x^2 + nx / sigma.sq.x[i]) 
      s1y = 1 / sqrt(1 / s0y^2 + ny / sigma.sq.y[i])
      
      m1x = m0x * s1x^2  / s0x^2  + Sx * s1x^2 / sigma.sq.x[i] 
      m1y= m0y * s1y^2  / s0y^2  + Sy * s1y^2 / sigma.sq.y[i]
      
      mu.x[i] = rnorm(1, m1x, s1x) 
      mu.y[i] = rnorm(1, m1y, s1y) 
    }
    
    res = data.frame(mu.x = tail(mu.x, nIter),
                     mu.y = tail(mu.y, nIter),
                     mu.diff = tail(mu.x - mu.y, nIter),
                     sigma.sq.x = tail(sigma.sq.x, nIter),
                     sigma.sq.y = tail(sigma.sq.y, nIter),
                     tstat = tail((mu.x - mu.y) / sqrt(sigma.sq.x / nx + sigma.sq.y / ny), nIter))
  }
  
  return(res)
}