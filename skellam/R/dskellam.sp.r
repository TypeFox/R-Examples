dskellam.sp <- function(x, lambda1, lambda2=lambda1, log=FALSE){
 # saddlepoint density (PMF) for Skellam distribution
  terms=1
    if (missing(x)|missing(lambda1)) stop("first 2 arguments are required")
    s <- log(0.5*(x+sqrt(x^2+4*lambda1*lambda2))/lambda1)# the saddlepoint
    K <- lambda1*(exp(s)-1)+lambda2*(exp(-s)-1)	# CGF(s)
    K2 <- lambda1*exp(s)+lambda2*exp(-s)        # CGF''(s)
    if (terms<1) {
        ret <- exp(K-x*s)/sqrt(2*pi*K2)         # saddlepoint density
    } else {
        c <- (1-((lambda1*exp(s)-lambda2*exp(-s))/K2)^2*5/3)/K2*0.125+1
        ret <- exp(K-x*s)/sqrt(2*pi*K2)*(1+c)*0.5
    }
    ret
}
