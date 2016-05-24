sunter <-
function (x, n) 
{
    N <- length(x)
    e = runif(N)
    x1 <- sort(x, decreasing = TRUE)
    if(n*x1[1]/sum(x1)>1) stop("There are some units with inclusion probability >1")
    t<-rev(cumsum(rev(x1)))
    kk0<-n*x1/t
    k0<-which(kk0>=1,arr.ind = FALSE)[1]
    kstar<-min(k0,N-n+1)
    nk <- numeric(N)
    I <- numeric(N)
    sam <- numeric(N)
    if (e[1] < n*x1[1]/t[1]) {
        I[1] <- 1
        sam[1] <- 1
    }
    k<-2
    nk[k]<-I[k-1]
    while (k<kstar & nk[k]<N-n+1) {
        d <- (n - nk[k])*x1[k]/t[k]
        if (e[k] <= d) {
            I[k] <- 1
            sam[k] <- cumsum(I[1:(k - 1)])[(k - 1)] + I[k]
        }
  k<-k+1
        nk[k] <- nk[k - 1] + I[k - 1]
    }
    while (nk[k]<n)
    {
    p<-(n-nk[k])/(N-k+1)
    if(e[k]<p)
    {I[k]<-1
     sam[k]<-cumsum(I[1:(k - 1)])[(k - 1)] + I[k]}
    k<-k+1
    nk[k]<-nk[k-1]+I[k-1]
    }
    samp <- rev(order(x))[which(sam != 0)]
    samp
}

