LM_stat <-
function(y,k)
{
    y <- as.matrix(y); y1 <- (y-mean(y))^2; n <- nrow(y); m <- mean(y)
    vr1 <- sum( (y-m)^2 )/n
  
     # use the filter function
    flt = filter(y, rep(1,k), method = "convolution")
    flt = flt[!is.na(flt)]
    summ = sum((flt - k * m)^2)

    vr2 <- summ/(n*k)
    vr <- vr2/vr1

    tem1 <- 2*(2*k-1)*(k-1)
    tem2 <- 3*k

    m1 <- sqrt(n)*(vr-1)/sqrt( tem1/tem2 )

    w <- 4*as.matrix((1-(1:(k-1))/k)^2,nrow=k-1)
    dvec <- matrix(NA, nrow=(k-1), ncol=1)
    for (j in 1:(k-1))
    { 
    dvec[j] <- sum(y1[(j+1):n] * y1[1:(n-j)])/( sum(y1)^2 )
    }
    summ <- crossprod(w,dvec)
    m2 <- sqrt(n)*(vr-1)*((n*summ)^(-.5) )
return(list(LM1=m1,LM2=m2))
}
