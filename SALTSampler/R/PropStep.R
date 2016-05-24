PropStep <- function(y, i, h) {
  #Given a logit-scaled simplex point y, this function 
  #draws a new logit-scaled simplex point.
  
  #Check input
  if (!exists("y")) {
    stop("y is not defined")
  }
  if (!is.numeric(y)){
    stop("y is not numeric")
  }
  
  if (!exists("i")) {
    stop("i is not defined")
  }
  if (!is.numeric(i)){
    stop("i is not numeric")
  }
  if (!is.vector(i)){
    stop("i is not a vector")
  }
  if (i%%1 != 0) {
    stop("i is not an integer")
  }
  if (length(i) != 1) {
    stop("i is not of length 1")
  }
  
  if (!exists("h")) {
    stop("h is not defined")
  }
  if (!is.numeric(h)){
    stop("h is not numeric")
  }
  if (!is.vector(h)){
    stop("h is not a vector")
  }
  
  #perturb ith logit
  ynew <- y
  ynew[i] <- y[i] + rnorm(1, 0, h)
  
  #Calculate logp and logq for old and new draws
  logpqOldNew <- LogPq(c(y[i], ynew[i]))
  
  #Log of the scaling value
  ls <- logpqOldNew$logq[2] - LogitSum(y[-i])
  
  #Logits of the rescaled simplex point
  ynew[-i] <- LogitScale(y[-i], ls)
  
  #Add detailed balance term
  dbt <- diff(logpqOldNew$logp) + (length(y) - 1)*diff(logpqOldNew$logq)
  attr(ynew, 'dbt') <- dbt
  
  #Return new logit-scaled point
  return(ynew)
}
