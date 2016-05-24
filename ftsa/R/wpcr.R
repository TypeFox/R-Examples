wpcr <- function(data, ncomp, method = c("ets", "rw"), h, beta = NULL, median = FALSE, transform = FALSE)
{
     n = ncol(data)
     q = matrix(,n,1)
     for(i in 1:n){
         q[i,] = beta * (1 - beta)^(i - 1)
     }
     weight = diag(rev(q))
     newdata = scale(t(data), scale = FALSE)
     if(median == TRUE){
        mdata = l1median(t(data))
     }
     if(median == FALSE){
        mdata = apply(data, 1, mean)
     }
     new = weight %*% newdata
     load = svd(new)$v[,1:ncomp]
     sco = newdata %*% load
     fore = matrix(NA, ncomp, h)
     if(method == "ets")
     {
         for(i in 1:ncomp)
         {
             fore[i,] = forecast(ets(sco[,i]), h = h)$mean
         }
     }
     if(method == "rw")
     {
         for(i in 1:ncomp)
         {
             fore[i,] = rwf(sco[,i], h = h, drift = FALSE)$mean
         }
     }
     if(transform == TRUE)
     {
         forecast = exp(load %*% fore + mdata)
     }
     if(transform == FALSE)
     {
         forecast = load %*% fore + mdata
     }
     return(forecast)
}