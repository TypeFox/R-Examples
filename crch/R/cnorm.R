## density
dcnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("dcnorm", x, mean, sd, left, right, log))
}

## distribution function
pcnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("pcnorm", q, mean, sd, left, right, lower.tail, log.p))
}

## random numbers
rcnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  rval <- rnorm(n) * sd + mean
  pmax(pmin(rval, right), left)
}

## quantiles
qcnorm <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  rval <- qnorm(p, lower.tail = lower.tail, log.p = log.p) * sd + mean
  pmax(pmin(rval, right), left)
}

## scores
scnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("scnorm_mu", x, mean, sd, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("scnorm_sigma", x, mean, sd, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
hcnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("hcnorm_mu", x, mean, sd, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("hcnorm_sigma", x, mean, sd, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("hcnorm_musigma", x, mean, sd, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}
