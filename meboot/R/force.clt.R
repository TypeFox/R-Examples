
force.clt <- function(x, ensemble)
{
  n <- nrow(ensemble)
  bigj <- ncol(ensemble)

  gm <- mean(x)  # desired grand mean
  #s <- sd(x)     # sd of original data
  s <- if (is.null(ncol(x)))
    sd(x) else apply(x, 2, sd)
  smean <- s/sqrt(bigj)  # desired standard deviation of means by CLT
  xbar <- apply(ensemble, 2, mean)
  sortxbar <- sort(xbar)  
  oo <- order(xbar)

  #now spread out the means as if they were from normal density
  #smallest mean should equal the normal quantile at prob=1/bigj+1
  #smallest second mean= the normal quantile at prob=2/bigj+1
  # . . .
  #LAST mean= the normal quantile at prob=bigj/(bigj+1)

  newbar <- gm + qnorm(1:bigj/(bigj+1))*smean

  #the above adjustement of means will change their sd
  # CLT says that their sd should be s / sqrt(n)
  #so we have recenter and rescale these revised means called newbar

  scn <- scale(newbar)  #first scale them to have zero mean and unit sd
  newm <- scn*smean+gm  #this forces the mean to be gm and sd=s / sqrt(n)

  meanfix <- as.numeric(newm - sortxbar)

  # we have lost the original order in sorting, need to go back 
  out <- ensemble
  for(i in 1:bigj)
    out[,oo[i]] <- ensemble[,oo[i]] + meanfix[i]

  if(any(is.ts(ensemble))){
    out <- ts(out, frequency=frequency(ensemble), start=start(ensemble))
    dimnames(out)[[2]] <- dimnames(ensemble)[[2]]
  }
  out
}
