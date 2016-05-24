WK_stat1 <-
function(y,k)
{
    y <- as.matrix(y)
    n <- nrow(y) 
    m <- mean(y)
    vr1 <- sum( (y-m)^2 )/n

    # use the filter function
    flt = filter(y, rep(1,k), method = "convolution")
    flt = flt[!is.na(flt)]
    summ = sum((flt - k * m)^2)

    vr2 <- summ/(n*k)
    vr <- vr2/vr1

    m1 <- sqrt(n)*(vr-1)
return(m1)
}
