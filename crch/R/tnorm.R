## density
dtnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("dtnorm", x, mean, sd, left, right, log))
}


## distribution function
ptnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("ptnorm", q, mean, sd, left, right, lower.tail, log.p))
}

## quantiles
qtnorm <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- pnorm((lower-mean)/sd, lower.tail = lower.tail) * (1 - p) + 
    p*pnorm((upper - mean)/sd, lower.tail = lower.tail)
  qnorm(p, lower.tail = lower.tail)*sd + mean
}

## random numbers
rtnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  qtnorm(runif(n), mean, sd, left = left, right = right)
}



## scores
stnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("stnorm_mu", x, mean, sd, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("stnorm_sigma", x, mean, sd, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("htnorm_mu", x, mean, sd, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("htnorm_sigma", x, mean, sd, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("htnorm_musigma", x, mean, sd, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}
