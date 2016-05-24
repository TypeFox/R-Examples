## ************************************************ ##
##       Logistic Split function                    ##
## ************************************************ ##

##  Initialization function

.logistic.init <- function(y, offset, params, wt)
{
  if(is.null(offset)) offset <- 0
  if(any(y != 1 && y != 0)) stop('Response most be 0/1')
  
  sfun <- function(yval, dev, wt, ylevel, digits)
  {
    paste("events = ", round(yval[,1]),
          ", coef = ", format(signif(yval[,2], digits)),
          ", deviance = ", format(signif(dev, digits)),
          sep=''
    )
  }
  tfun <- function (yval, dev, wt, ylevel, digits, n, use.n) {
    paste(format(signif(yval[,1]/n,digits)),"\n",
          format(signif(yval[,2], digits)),
          sep = '')}
  
  
  environment(sfun) <- .GlobalEnv
  environment(tfun) <- .GlobalEnv
  
  list(y = cbind(y, offset), params = 0, numresp = 2, numy = 2, summary = sfun, text=tfun)
}

##	Evaluation function

.logistic.eval <- function(y, wt, params)
{
  tfit <- glm(y[,1]~offset(y[,2]), binomial, weights = wt)
  list(label = c(sum(y[,1]), tfit$coef), deviance = tfit$deviance)
}

##	Spliting function

.logistic.split <- function(y, wt, x, params, continuous)
{
  if(continuous)
  {
    n <- nrow(y)
    goodness <- double(n-1)
    direction <- goodness
    temp <- rep(0, n)
    for(i in 1:(n-1))
    {
      temp[i] <- 1
      if(x[i] != x[i+1])
      {
        tfit <- glm(y[,1]~temp + offset(y[,2]), binomial, weights = wt)
        goodness[i] <- tfit$null.deviance - tfit$deviance
        direction[i] <- sign(tfit$coef[2])
      }
    }
  }
  else
  {
    x = x[,drop=TRUE]
    tfit <- glm(y[,1]~factor(x)+offset(y[,2])-1, binomial, weights = wt)
    ngrp <- length(tfit$coef)
    direction <- rank(rank(tfit$coef) + runif(ngrp, 0, 0.1)) 
    
    xx <- direction[match(x, sort(unique(x)))] 
    goodness <- double(length(direction) - 1)
    for(i in 1:length(goodness))
    {
      tfit <- glm(y[,1]~I(xx > i) + offset(y[,2]), binomial, weights = wt)
      goodness[i] <- tfit$null.deviance - tfit$deviance
    }
  }
  
  list(goodness = goodness, direction = direction)
}
