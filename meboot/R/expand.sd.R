
expand.sd <- function(x, ensemble, fiv=5)
{
  sdx <- if (is.null(ncol(x)))
    sd(x) else apply(x, 2, sd)

  sdf <- c(sdx, apply(ensemble, 2, sd))
  sdfa <- sdf/sdf[1]  # ratio of actual sd to that of original data
  sdfd <- sdf[1]/sdf  # ratio of desired sd to actual sd

  # expansion is needed since some of these are <1 due to attenuation
  mx <- 1+(fiv/100)
  # following are expansion factors
  id <- which(sdfa < 1)
  if (length(id) > 0)
    sdfa[id] <- runif(n=length(id), min=1, max=mx)

  # multiply each column by expansion factor
  #
  # June 20, 2014
  # Do not apply if sdfd[-1]*sdfa[-1] < 1 since 
  # diag(sdfd[-1]*sdfa[-1]) would be empty
  # Thanks to Benjamin Tyner for reporting this issue
  #
  # July 26, 2014 thanks to Régis Hank for reporting a bug
  # "id" is necessary

  sdfdXsdfa <- sdfd[-1]*sdfa[-1]
  id <- which(floor(sdfdXsdfa) > 0)
  if (length(id) > 0) {
    ensemble[,id] <- ensemble[,id] %*% diag(sdfdXsdfa[id])
  } else
  
  # changed June 20, 2014
  #if(any(is.ts(ensemble))){
  if(is.ts(x)){
    ensemble <- ts(ensemble, frequency=frequency(x), start=start(x))
    #dimnames(out)[[2]] <- dimnames(ensemble)[[2]]
  }  
  ensemble
}
