`TacvfMA` <-
function(theta, lag.max = 20)
{
    if(length(theta) == 0)
        if(lag.max >= 0) 
            return(c(1, numeric(lag.max+1)))
        else 
            stop("maxlag invalid")
    maxlagp1 <- lag.max+1
    g<-rep(0, maxlagp1)
    th<-c(-1,theta)
    qth <- length(th)
    x<-c(th, rep(0,qth))
    A <- matrix(0, qth, qth)
    B<-matrix(x[abs(col(A) + row(A)) - 1], qth, qth)
    g1<-c(B%*%th)
    if (length(g1)<maxlagp1)
        g[1:qth]<-g1
    else
        g<-g1[1:maxlagp1]
    g
}

