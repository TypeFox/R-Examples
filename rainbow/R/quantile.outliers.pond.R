`quantile.outliers.pond` <- function(data, dfunc = depth.RP, nb = 200, suav = 0.05, ...)
{
    functions = t(data$y)
    n <- dim(functions)[1]
    m <- dim(functions)[2]
    if(is.null(n) && is.null(m)) 
       stop("I do not have a matrix.")
       d = dfunc(data,...)$prof
       quantiles <- numeric(nb)
       vv = var(functions)
       pr = d / sum(d)
       for(i in 1:nb){
           bsample <- functions[sample(1:n, size = n, replace = T, prob = pr),]
           if(suav>0){
              bsample <- bsample + mvrnorm(n = n, rep(0, m), vv * suav)
           }
           bsample = fts(1:dim(bsample)[1], bsample)
           d = dfunc(bsample, ...)$prof
           quantiles[i] <- quantile(d, probs = 0.01, type = 8)
       }
     return(quantiles)
}

