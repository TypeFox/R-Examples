#' @title Likelihood Functions for Generating OpenBUGS Model File
#' @aliases bugsLikelihood
#' @description Likelihood functions for generating OpenBUGS model file.
#' 
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' @return A string to be written in the BUGS model file.
#' 
#' @export
#' @seealso \link{qrBugs}
#' @examples
#' bugsLikelihood('t2','t2')

bugsLikelihood <- function(fd, sd) {
  # logit-XX
  if (fd == "logit") {
    
    if (sd == "logistic") 
      loglik <- "logLike[i] <-  logpart1[i] + logpart2[i] \n    logpart1[i] <- ((mu[i] - sigma[i] * log(sigma[i])) - (-1 + sigma[i]) * log(-1 + 1/y[i]))/sigma[i]  \n    logpart2[i] <- -(2 * sigma[i] * log(y[i] +exp(mu[i]/sigma[i]) * pow((-1 + 1/y[i]),\n    (1/sigma[i])) * y[i]))/sigma[i]"
    
    if (sd == "cauchy") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] + logpart3[i] \n    logpart1[i] <-  (mu[i] + sigma[i] * log(pi) - sigma[i] * log(sigma[i]) +(1/tan(pi * y[i])) )/sigma[i]  \n    logpart2[i] <- -(2 * sigma[i] * log(1 + exp((mu[i] + (1/tan(pi * y[i])))/sigma[i])))/sigma[i]\n    logpart3[i] <- (2 * sigma[i] * log(1/sin(pi * y[i])))/sigma[i]"
    
    if (sd == "t2") 
      loglik <- "case1[i] <- (step(y[i] - 1/2) - 1/2)*2 #case1[i] = 1 if y[i]>= 1/2; case1[i] = -1 if y[i]< 1/2\n    a1[i] <- case1[i]*sqrt(pow(1 - 2 * y[i],2)/(2 * (1 - y[i]) * y[i]))  \n    case2[i] <- step(y[i] - 1/2)*step(1/2 - y[i]) #case2[i] = 1 if y[i] = 1/2; case2[i] = 0 if y[i] != 1/2    \n    a2[i] <- case2[i]*2 * sqrt(2) + (1-case2[i])*case1[i]*(2 * y[i] - 1)/sqrt(8 * pow(1 - 2 * y[i],2) * pow((1 - y[i]) * y[i],3))\n    logLike[i]  <- log((exp((mu[i] + a1[i])/sigma[i]) * a2[i])/(pow(exp(mu[i]/sigma[i]) + exp(a1[i]/sigma[i]),2) * sigma[i]))"
    
    
    if (sd == "burr8") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] + logpart3[i] \n    logpart1[i] <-  log(1/sin(pi * y[i]))\n    logpart2[i] <- ((mu[i] + sigma[i] * log(pi) - sigma[i] * log(sigma[i])) + log(tan((pi * y[i])/2)))/sigma[i]    \n    logpart3[i] <- -(2 * sigma[i] * log(exp(mu[i]/sigma[i]) + pow(tan((pi * y[i])/2),1/sigma[i])))/sigma[i]"
    
  }
  
  # archsinh-XX
  if (fd == "arcsinh") {
    
    if (sd == "cauchy") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] + logpart3[i] \n    logpart1[i] <-  log(pi) + arcsinh((mu[i] + (1/tan(pi * y[i])))/sigma[i])  \n    logpart2[i] <- -2 * log(1 + exp(arcsinh((mu[i] + (1/tan(pi * y[i])))/sigma[i])))+ 2 * log(1/sin(pi * y[i]))\n    logpart3[i] <- -(1/2) * log(pow(mu[i], 2) + pow(sigma[i], 2) + 2 * mu[i] * 1/tan(pi * y[i]) +pow(1/tan(pi * y[i]),2))"
    
    
    
    if (sd == "t2") 
      loglik <- "\n    case1[i] <- (step(y[i] - 1/2) - 1/2)*2 #case1[i] = 1 if y[i]>= 1/2; case1[i] = -1 if y[i]< 1/2\n    a1[i] <- case1[i]*sqrt(pow((1 - 2 * y[i]),2)/(2 * (1 - y[i]) * y[i]))  \n    \n    case2[i] <- step(y[i] - 1/2)*step(1/2 - y[i]) #case2[i] = 1 if y[i] = 1/2; case2[i] = 0 if y[i] != 1/2   \n    a2[i] <- case2[i]*2 * sqrt(2) + (1-case2[i])*case1[i]*(2 * y[i] - 1)/sqrt(8 * pow((1 - 2 * y[i]),2) * pow(((1 - y[i])* y[i]),3))\n    \n    logpart1[i] <- exp(arcsinh((-mu[i] + a1[i])/sigma[i])) * a2[i]\n    logpart2[i] <- pow((1 + exp(arcsinh((-mu[i] + a1[i])/sigma[i]))),2) * \n    sigma[i] *  sqrt(1 + pow((mu[i] - a1[i]), 2)/pow(sigma[i], 2))\n    logLike[i]  <- log(logpart1[i]/logpart2[i])\n    "
    
    if (sd == "burr8") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] + logpart3[i] \n    logpart1[i] <- log(pi) - log(sigma[i]) + arcsinh((log(tan((pi * y[i])/2))-mu[i])/sigma[i])\n    logpart2[i] <- -2 * log(1 +  exp(arcsinh((-mu[i] + log(tan((pi * y[i])/2)))/sigma[i])))    \n    logpart3[i] <- log(1/sin(pi * y[i])) - (1/2) * log(1 + pow((mu[i] - log(tan((pi * y[i])/2)))/sigma[i], 2))"
    
  }
  
  # t2-XX
  if (fd == "t2") {
    
    if (sd == "cauchy") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] + logpart3[i] \n    logpart1[i] <-  log(pi) + log(sigma[i])   \n    logpart2[i] <- -log(pow(mu[i], 2) + 2 * pow(sigma[i], 2) + (2 * mu[i]/tan(pi * y[i])) + pow((1/tan(pi * y[i])),2)) \n    logpart3[i] <- -(1/2) * log(2 + pow((mu[i] +( 1/tan(pi * y[i]))),2)/pow(sigma[i], 2)) + 2 * log(1/sin(pi * y[i]))"
    
    
    
    if (sd == "t2") 
      loglik <- "\n    logLike[i] <-  2 * log(sigma[i])  -(3/2) * log(logpart1[i] + logpart2[i] + logpart3[i])\n    logpart1[i] <-  1 + 2 * sqrt(2) * mu[i] * sqrt(1 - y[i]) * sqrt(y[i])   \n    logpart2[i] <-  2 * (-2 + pow(mu[i], 2) + 2 * pow(sigma[i], 2)) * y[i]\n    logpart3[i] <- -4 * mu[i] * sqrt(2 - 2 * y[i]) * pow(y[i], 1.5) - 2 * (-2 + pow(mu[i], 2) + 2 * pow(sigma[i], 2)) * pow(y[i], 2)"
    
    if (sd == "burr7") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] \n    logpart1[i] <- -log(2) + 2 * log(sigma[i])- log(1 - y[i]) - log(y[i])\n    logpart2[i] <-  -(3/2) * log(pow(mu[i], 2) + 2 * pow(sigma[i], 2) + 2 * mu[i] * arctanh(1 - 2 * y[i]) + \n    pow(arctanh(1 - 2 * y[i]),2))"
    
    if (sd == "burr8") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i]\n    logpart1[i] <- log(pi) + 2 * log(sigma[i]) + log(1/sin(pi * y[i]))\n    logpart2[i] <- -(3/2) * log(pow(mu[i], 2) + 2 * pow(sigma[i], 2) - 2 * mu[i] * log(tan((pi * y[i])/2)) + pow(log(tan((pi * y[i])/2)), 2))"
    
  }
  
  # burr8-XX
  if (fd == "burr8") {
    
    if (sd == "cauchy") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] \n    logpart1[i] <-  (mu[i] + sigma[i] * log(2) - sigma[i] * log(sigma[i]) + 1/tan(pi * y[i]))/sigma[i]  \n    logpart2[i] <- -log(1 + exp((2 * (mu[i] + 1/tan(pi * y[i])))/sigma[i]))+2 * log(1/sin(pi * y[i]))"
    
    
    if (sd == "t2") 
      loglik <- "\n    logLike[i]  <- log(logpart1[i]/logpart2[i])\n    logpart1[i] <- 1/cosh(((1 - 2 * y[i])/(sqrt(2) * sqrt((1-y[i]) * y[i])) + mu[i])/sigma[i])\n    logpart2[i] <- 2 * sqrt(2) * pi * pow(((1- y[i])*y[i]), 3/2) * sigma[i]\n    "
    
    
    if (sd == "burr7") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i] \n    logpart1[i] <- -log(2 * pi * y[i] * sigma[i] - 2 * pi * pow(y[i], 2) * sigma[i])\n    logpart2[i] <-  log(1/cosh((mu[i] + arctanh(1 - 2 *  y[i]))/sigma[i]))"
    
    
    if (sd == "burr8") 
      loglik <- "\n    logLike[i] <-  logpart1[i] + logpart2[i]\n    logpart1[i] <-  (mu[i] + sigma[i] * log(2) + sigma[i] * log(1/sin(pi * y[i])) + log(tan(pi * y[i]/2)))/sigma[i] \n    logpart2[i] <- -log(exp((2 * mu[i])/sigma[i]) * sigma[i] + sigma[i] * pow(tan((pi * y[i])/2), 2/sigma[i]))"
    
  }
  
  if (fd == "km" | sd == "km") {
    loglik <- "\n    logLike[i] <- log(mu[i]) + log(sigma[i]) + (mu[i] - 1)*log(y[i]) + 
    (sigma[i] - 1)*log(1 - pow(y[i],mu[i]))"
  }
  
  
  loglik
} 
