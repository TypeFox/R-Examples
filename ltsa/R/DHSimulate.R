"DHSimulate" <-
function(n, r, ReportTestOnly=FALSE, rand.gen=rnorm, ...){
    m <- length(r)
    N <- 2^ceiling(log(n-1,base=2))
    if (n == m-1)
            acvf <- r
    else
            acvf <- c(r,rep(0,N-m+1))
    d <- c(acvf,rev(acvf[-1])[-1])
    g <- Re(fft(d))
    test <- any(g<0)
    if (ReportTestOnly)
         return(!test)
    if (test)
        stop("Davies-Harte nonnegativity condition not valid")
    Z<-complex(real = rand.gen(N-1, ...), imaginary = rand.gen(N-1, ...))
    Z2<-2+sqrt(2)*rand.gen(2, ...)
    Z<-c(Z2[1],Z,Z2[2],Conj(rev(Z)))
    X<-Re(fft(sqrt(g)*Z,inverse=TRUE))/sqrt(2*N)
    z<-X[1:n]/sqrt(2)
    z
}
