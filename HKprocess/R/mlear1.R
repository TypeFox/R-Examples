mlear1 <- function(data,interval = c(-0.9999,0.9999),
                   tol = .Machine$double.eps^0.25) {
    ar1likelihood <- function(f1,x) {
        if (f1 <= -1) return(-Inf)
        if (f1 >= 1) return(-Inf)
        nx <- as.integer(length(x))
        r <- f1^(0:(nx - 1))
        EPS<-.Machine$double.eps

        # Call the likelihoodfunction.c from the C library
        out<-.C("likelihoodfunction",as.double(r),nx,as.double(x),nx,EPS,
        tr = array(0,dim = c(1,4)),fault = as.integer(1),PACKAGE = "HKprocess")

        return(c((out$tr)[1],(out$tr)[2],(out$tr)[3],(out$tr)[4]))
    }
    
    # Set the function f which will be optimized to obtain the estimate of H
    f <- function(f1,x) (ar1likelihood(f1,x))[1]

    y1 <-as.double((optimize(f,interval = interval,x = data,maximum = TRUE,
     tol = tol))[1]) # y1 is the estimate of H
    
    y2 <- as.vector((ar1likelihood(y1,data))[3:4])
    # y2 is a vector with the estimates of mu and sigma
    return(setNames(c(y2,y1),c("mu_estimate","sigma_estimate","phi_estimate")))
}