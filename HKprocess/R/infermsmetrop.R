infermsmetrop <- function(Hbayes,data) {
    # Some manipulations to get the scale of the inv-gamma distribution for the
    # simulation of sigma^2, using eq. 9 in Tyralis and Koutsoyiannis (2014) 
    logpHx <- function(H,x) {
        if (H <= 0) return(-Inf)
        if (H >= 1) return(-Inf)
        nx <- as.integer(length(x))
        maxlag <- nx - 1
        r <- acfHKp(H,maxlag)
        EPS<-.Machine$double.eps

        out <- .C("logpHx",as.double(r),nx,as.double(x),nx,EPS,
        tr = array(0,dim=c(1,5)),fault = as.integer(1),PACKAGE = "HKprocess")

        return(c((out$tr)[1],(out$tr)[2],(out$tr)[3],(out$tr)[4],(out$tr)[5]))
    }
    g <- function(H,x = data) (logpHx(H,data))
    w <- as.vector(sapply(Hbayes,g)) # w contains the required scales
    size = length(data)
    sizeH = length(Hbayes)
    shape = 0.5 * (size - 1)
    minfer = vector()
    sinfer2 = vector()
    # Simulation of eqs. 8,9 in Tyralis and Koutsoyiannis (2014)
    for (i in 1:sizeH) {
        sinfer2[i] = rinvgamma(1,shape,scale = w[5*i])
        minfer[i] = rnorm(1,mean = w[5*i - 2],sd = sqrt(sinfer2[i]*w[5*i - 1]))}
    
    return(matrix(c(minfer,sinfer2),nrow = sizeH,ncol = 2,dimnames = list(NULL,
    c("mu_sample","sigma_sq_sample"))))
}