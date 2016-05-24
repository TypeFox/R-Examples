VR.minus.1 <-
function(y,kvec)
{
    coe <- AR1(y)$ALPHA
    T <- length(y)
    lq <- ABEL1Q(T,coe)    
    vrsum <- 1
    for (i in 1:(T-1))
    {
        sum1 <- 0
        for (t in 1:(T-i))
        {
        sum1 <- sum1 + y[t]*y[t+i]
        }
    sum1 <- sum1/(sum(y^2))
    vrsum <- vrsum + 2*kfunc(i/lq)*sum1
    }
vr.auto <- (vrsum - 1)

    y <- as.matrix(y)
    n <- nrow(y) 
    m <- mean(y)
    vr1 <- sum( (y-m)^2 )/n
    
    mq <- numeric()
    for (i in 1:length(kvec))
    {
    k <- kvec[i]
    
    # use the filter function
    flt = filter(y, rep(1,k), method = "convolution")
    flt = flt[!is.na(flt)]
    summ = sum((flt - k * m)^2)

    vr2 <- summ/(n*k)
    vr <- vr2/vr1
    
    mq <- c(mq,(vr-1))
    }  
   # rownames(mq) <- paste("k",kvec,sep=""); colnames(mq) <- "|VR-1|"
return(list(VR.auto=vr.auto,Holding.Periods=kvec,VR.kvec=mq))
}
