mleHK <- function(data,interval = c(0.0001,0.9999),
                  tol = .Machine$double.eps^0.25) {
    HKlikelihood <- function(H,x) {
        if (H <= 0) return(-Inf)
        if (H >= 1) return(-Inf)
        nx <- as.integer(length(x))
        maxlag <- nx - 1
        r <- acfHKp(H,maxlag)
        EPS<-.Machine$double.eps

        # Call the likelihoodfunction.c from the C library
        out<-.C("likelihoodfunction",as.double(r),nx,as.double(x),nx,EPS,tr =
        array(0,dim = c(1,4)),fault = as.integer(1),PACKAGE = "HKprocess")
        return(c((out$tr)[1],(out$tr)[2],(out$tr)[3],(out$tr)[4]))
    }
    # Set the function f which will be optimized to obtain the estimate of H
    f <- function(H,x) (HKlikelihood(H,x))[1]

    y1 <-as.double((optimize(f,interval = interval,x = data,maximum = TRUE,
    tol = tol))[1])
    # y1 is the estimate of H

    y2<-as.vector((HKlikelihood(y1,data))[3:4])
    # y2 is a vector with the estimates of mu and sigma
    return(setNames(c(y2,y1),c("mu_estimate","sigma_estimate","H_estimate")))
}