## density
dclogis <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("dclogis", x, mean, sd, left, right, log))
}

## distribution function
pclogis <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  with(input, .Call("pclogis", q, mean, sd, left, right, lower.tail, log.p))
}

## random numbers
rclogis <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  rval <- rlogis(n) * sd + mean
  pmax(pmin(rval, right), left)
}

## quantiles
qclogis <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  rval <- qlogis(p, lower.tail = lower.tail, log.p = log.p) * sd + mean
  pmax(pmin(rval, right), left)
}

## scores
sclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("sclogis_mu", x, mean, sd, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("sclogis_sigma", x, mean, sd, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

# Hessian
hclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("hclogis_mu", x, mean, sd, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("hclogis_sigma", x, mean, sd, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("hclogis_musigma", x, mean, sd, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

### density
#dclogis <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
#  x <- data.frame(x = x, mean, sd)$x
#  ifelse(x <= left, plogis((left - mean)/sd, log.p = log), 
#  ifelse(x >= right, plogis((right - mean)/sd, log.p = log, lower.tail = FALSE), 
#  dlogis((x - mean)/sd, log = log)/sd^(1 - log) - log(sd) * log))
#}

### distribution function
#pclogis <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
#  left = -Inf, right = Inf) {
#  q <- data.frame(q = q, mean, sd)$q
#  if(lower.tail){ 
#    ifelse(q < left, 0, 
#    ifelse(q >= right, 1, 
#    plogis((q-mean)/sd, log.p = log.p)))
#  } else {
#    ifelse(q <= left, 1, 
#    ifelse(q > right, 0, 
#    plogis((q-mean)/sd, lower.tail = FALSE, log.p = log.p)))
#  }
#}


### scores
#sclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf) {
#  if(!is.character(which))
#    which <- c("mu", "sigma")[as.integer(which)]
#  which <- tolower(which)
#  score <- NULL
#  x <- data.frame(x = x, mean, sd)$x
#  dxm <- x - mean
#  dlm <- left - mean
#  drm <- right - mean 
#  sd2 <- sd^2
#  millsl <- dlogis(dlm/sd)/sd/plogis(dlm/sd)
#  millsr <- dlogis(drm/sd)/sd/plogis(drm/sd, lower.tail = FALSE)
#  for(w in which) {
#    if(w == "mu")
#      score2 <- ifelse(x <= left, - millsl,
#        ifelse(x >= right, millsr,
#        (1 - 2 * plogis(-dxm/sd))/sd))
#    if(w == "sigma")
#      score2 <- ifelse(x <= left, - millsl * dlm/sd,
#        ifelse(x >= right, millsr * drm/sd,
#        (1 - 2 * plogis(-dxm/sd))*dxm/sd2 - 1/sd))
#    score <- cbind(score, score2)
#  }
#  if(is.null(dim(score)))
#    score <- matrix(score, ncol = 1)
#  colnames(score) <- paste("d", which, sep = "")
#  score
#}

### Hessian
#hclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf)
#{
#  if(!is.character(which))
#    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
#  which <- tolower(which)
#  hess <- list()
#  x <- data.frame(x = x, mean, sd)$x
#  sd2 <- sd^2
#  dlm <- left - mean
#  drm <- right - mean
#  dxm <- x - mean
#  millsl <- dlogis(dlm/sd)/sd/plogis(dlm/sd)
#  millsr <- dlogis(drm/sd)/sd/plogis(drm/sd, lower.tail = FALSE)
#  scorel <- sclogis(left, mean, sd, which = "mu", left = -Inf, right = Inf)
#  scorer <- sclogis(right, mean, sd, which = "mu", left = -Inf, right = Inf)
#  score  <- sclogis(x, mean, sd, left = -Inf, right = Inf)
#  for(w in which) {
#    if(w == "mu")
#      hess[[w]] <- 
#          ifelse(  x <= left, - scorel * millsl - millsl^2,
#            ifelse(x >= right,  scorer * millsr - millsr^2,
#                                - 2/sd2 * dlogis(dxm/sd)))
#    if(w == "sigma")
#      hess[[w]] <- 
#        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^2/sd2*scorel)*millsl-
#          millsl^2*dlm^2/sd2,
#          ifelse(x >= right, (- 2 * drm/sd2 + drm^2/sd2*scorer)*millsr-
#            millsr^2*drm^2/sd2,
#                             - score[,"dmu"]*dxm/sd2 - 2 * dxm^2/sd2^2 * 
#                              dlogis(dxm/sd) - score[, "dsigma"]/sd))

#    if(w %in% c("mu.sigma", "sigma.mu"))## not correct
#      hess[[w]] <- 
#        ifelse(x <= left,    (  1/sd - dlm/sd*scorel) * millsl - 
#          dlm/sd * millsl^2,
#          ifelse(x >= right, (- 1/sd + drm/sd*scorer) * millsr - 
#            drm/sd * millsr^2,
#                             -score[, "dmu"]/sd - 2*dxm/sd^3*dlogis(dxm/sd)))
#  }

#  hess <- do.call("cbind", hess)
#  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
#  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
#  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
#  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
#  hess
#}
