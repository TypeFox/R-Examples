## density
dtlogis <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("dtlogis", x, mean, sd, left, right, log))
}


## distribution function
ptlogis <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("ptlogis", q, mean, sd, left, right, lower.tail, log.p))
}

## quantiles
qtlogis <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- plogis((lower-mean)/sd, lower.tail = lower.tail) * (1 - p) + 
    p*plogis((upper - mean)/sd, lower.tail = lower.tail)
  qlogis(p, lower.tail = lower.tail)*sd + mean
}

## random numbers
rtlogis <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  qtlogis(runif(n), mean, sd, left = left, right = right)
}

## scores
stlogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("stlogis_mu", x, mean, sd, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("stlogis_sigma", x, mean, sd, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htlogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("htlogis_mu", x, mean, sd, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("htlogis_sigma", x, mean, sd, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("htlogis_musigma", x, mean, sd, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}
