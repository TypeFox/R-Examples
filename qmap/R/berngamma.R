dberngamma <- function(x, prob, scale, shape){
  if(length(prob)==1) prob <- rep(prob, length(x))
  if(length(scale)==1) scale <- rep(scale, length(x))
  if(length(shape)==1) shape <- rep(shape, length(x))
  d <- 1-prob
  d[x>0] <- prob[x>0]*dgamma(x[x>0], scale=scale[x>0], shape=shape[x>0])
  d
}


pberngamma <-  function(q, prob, scale, shape){
  if(length(prob)==1) prob <- rep(prob, length(q))
  if(length(scale)==1) scale <- rep(scale, length(q))
  if(length(shape)==1) shape <- rep(shape, length(q))
  p <- 1-prob
  p[q>0] <- 1-prob[q>0]+prob[q>0]*pgamma(q[q>0], scale=scale[q>0],
                                         shape=shape[q>0])
  p 
}


qberngamma <- function(p, prob, scale, shape){
  if(length(prob)==1) prob <- rep(prob, length(p))
  if(length(scale)==1) scale <- rep(scale, length(p))
  if(length(shape)==1) shape <- rep(shape, length(p))
  q <- rep(0, length(p))
  cases <- p > (1-prob)
  q[cases] <- qgamma((prob[cases]+p[cases]-1)/prob[cases],
                     scale=scale[cases], shape=shape[cases])
  q
}


rberngamma <- function(n, prob, scale, shape){
  if(max(length(prob), length(scale), length(shape)) > 1)
    stop("parameters must be of length 1")
  p <- runif(n)
  q <- rep(0, length(p))
  cases <- p > (1-prob)
  q[cases] <- rgamma(sum(cases), scale=scale, shape=shape)
  q
}
