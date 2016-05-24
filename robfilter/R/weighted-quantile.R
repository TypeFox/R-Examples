weighted.quantile.top_neu = function(dat,weights,p){
    ord <- order(dat)
    weights <- weights[ord]
    n <- length(weights)
    dat <- dat[ord]
    q <- sum(weights) * p
    ind <- n+1-which(cumsum(rev(weights))>=q)[1]
    dat[ind]
}
weighted.quantile.top <- function(dat, weights, p){
  weights <- weights[order(dat)]
  n       <- length(weights)
  dat     <- sort(dat)
  q       <- sum(weights)*p
  m<- w   <-0
  for (i in n : 1){
    w <-  w+weights[i]
    if (w >= q) {m<- dat[i] ; break}
  }
  m
 }

dat <- rnorm(20)
weights <- rnorm(20)
p <- runif(1,0,1)
weighted.quantile.top_neu(dat,weights,p)
weighted.quantile.top(dat,weights,p)

weighted.quantile_neu = function(dat,weights,p){
    ord <- order(dat)
    weights <- weights[ord]
    n <- length(weights)
    dat <- dat[ord]
    q <- sum(weights) * p
    ind <- which(cumsum(weights)>=q)[1]
    dat[ind]
}
weighted.quantile <- function(dat, weights, p){
  weights <- weights[order(dat)]
  dat     <- sort(dat)
  q       <- sum(weights)*p
  m <- w  <-0
  for (i in 1 : length(weights)){
    w <- w+weights[i]
    if (w >= q) {m<- dat[i] ; break}
  }
  m
}

dat <- rnorm(20)
weights <- rnorm(20)
p <- runif(1,0,1)
weighted.quantile_neu(dat,weights,p)
weighted.quantile(dat,weights,p)
