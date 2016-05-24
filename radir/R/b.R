b <-
function(x, pos=1, envir=as.environment(pos))
{ 
  u2 <- get("u.local2", envir=envir)(x)
  v2 <- get("v.local2", envir=envir)[[1]](x)(x)
  pr <- c(attr(v2, "gradient"))
  return(u2/(pr %*% cov %*% as.matrix(pr)))
}
