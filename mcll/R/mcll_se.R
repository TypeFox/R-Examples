# Standard error estimation

mcll_se <- function(data, par,  H.prior , alp=0.7,        
        method="Nelder-Mead",  lower = -Inf, upper = Inf, control=list() ) {

    # data: posterior samples 
    # par: parameter estimates 
    # H.prior: Hessian matrix of log prior 
    # alp: alpha smoothing parameter  
    # method: optimization method 
    # lower, upper: boundary values  for the "L-BFGS-B" method
    # control: control options for optim 

    # covariance matrix    
    Sig <- cov(data)    
    
    # lower triangle 
    ch <- t(chol(Sig))
    
    # mean 
    bb <- colMeans(data) 
    
    # orthonormalization 
    trans <- function(vec, ch, bb) {
        solve(ch)%*% t(matrix( as.matrix(vec), nrow=1, ncol=length(bb)) - matrix(bb, nrow=1, ncol=length(bb) ) )  
    }
    
    data.t <- t(apply(data, 1, trans, ch=ch, bb=bb))

    par.t <- as.vector(t(solve(ch)%*%matrix(as.numeric(par) - bb, length(par),1)))

    # quadrature locations and weights 
    out <- gauss.quad(30,"legendre",-1,1)
    nodes=out$nodes ; weights=out$weights
    
    # optimization options 
    npar <- length(par.t)  # dimension 
    start <- c( (-0.5*npar*log(2*pi)), rep(0,npar), rep(-1,npar))
    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
        length(start) ), ndeps = rep.int(0.001, length(start)), maxit = 100L, abstol = -Inf, 
        reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
        gamma = 2, REPORT = 10, type = 1, lmm = 5, factr = 1e+07, 
        pgtol = 0, tmax = 10, temp = 10)    

    con[(namc <- names(control))] <- control

    try(result <- mcll_coeff(start=start, data.t=data.t,  par.t= par.t, alp=alp,  
                    nodes=nodes, weights=weights, 
                    method=method, lower = lower, upper = upper,  control=con )   ) 
    
    # return values of mcll_coeff 
    pol_coeff <- as.vector(result$a_coeff)
    convergence <- result$convergence
    value <- result$value
    counts <- result$counts
    message <- result$message 
    
    if (convergence !=0)
     warning("Convergence might not have been reached.")
    
    # second derivative of log posterior 
    a_coeff <- pol_coeff[(length(pol_coeff)-(npar-1)):length(pol_coeff)]
    der2 <- as.numeric(a_coeff)
    log.post.D <- diag(der2)    
        
    ## second derivative of the log-likelihood 
        
    lik.D <-log.post.D - H.prior 
        
    # standard error 
    cov.t <- ch %*% solve(-lik.D) %*% t(ch)
    se <- sqrt(diag(cov.t))
        
    return(se) 

}
