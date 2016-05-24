`quantile.outliers.trim` <- function(data, dfunc = depth.mode, trim = 0.1, nb = 200, suav = 0.05,...)
{
  functions = t(data$y)
  n <- dim(functions)[1]
  m <- dim(functions)[2]
  if(is.null(n) && is.null(m)) 
     stop("I do not have a matrix")
  d = dfunc(data, trim = trim,...)$prof
  rid <- rank(d, ties.method = "first")
  num.boot <- floor(trim * n)
  sample.trim <- functions[rid >= num.boot,]
  cuantiles <- numeric(nb)
  vv = var(functions)
  for(i in 1:nb){
      bsample <- sample.trim[sample(1:(n - num.boot), size = n, replace = T),]
      if(suav>0){
         bsample <- bsample + mvrnorm(n = n, rep(0,m), vv * suav)
      }
      bsample = fts(1:dim(bsample)[1], bsample)
      d = dfunc(bsample,...)$prof
      cuantiles[i] <- quantile(d, probs = 0.01, type = 8)
  }
  return(cuantiles)
}

