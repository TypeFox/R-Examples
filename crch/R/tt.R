## density
dtt <- function(x, mean = 0, sd = 1, df, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("dtt", x, mean, sd, df, left, right, log))
}


## distribution function
ptt <- function(q, mean = 0, sd = 1, df, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("ptt", q, mean, sd, df, left, right, lower.tail, log.p))
}

## quantiles
qtt <- function(p, mean = 0, sd = 1, df, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  p <- pt((lower-mean)/sd, df = df) * (1 - p) + p*pt((right - mean)/sd, df = df)
  qt((p - mean)/sd, df = df, lower.tail = lower.tail)
}

## random numbers
rtt <- function(n, mean = 0, sd = 1, df, left = -Inf, right = Inf) {
  qtt((runif(n) - mean)/sd, df = df, left = left, right = right)
}

## scores
stt <- function(x, mean = 0, sd = 1, df, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("stt_mu", x, mean, sd, df, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("stt_sigma", x, mean, sd, df, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htt <- function(x, mean = 0, sd = 1, df, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("htt_mu", x, mean, sd, df, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("htt_sigma", x, mean, sd, df, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("htt_musigma", x, mean, sd, df, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}
