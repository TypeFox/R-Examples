# Parameter estimation 

lik_manual <- function(x,  data.t, ch, bb, prior.func, alp,
            method2,  lower2 , upper2 , control2 ) {

    # data.t: transformed data (of the posterior samples of the parameters)
    # ch: lower triangular matrix of the cholesky for the original posterior samples 
    # bb: means of the original posterior samples 
    # prior.func: prior function (user specified) that provies a log prior density 
    # alp: alpha smoothing parameter 
          
    vec.t <- as.numeric(x)  # transformed scale 
    npar <- length(vec.t) # dimension 

    # transfomation back to the original scale 
    vec <- t((ch%*%vec.t) +bb)
    
    # quadrature locations and weights 
    out <- gauss.quad(30,"legendre",-1,1)
    nodes=out$nodes ; weights=out$weights
    
    # starting values 
    npar <- length(vec)  # dimension 
    aa.start <- c( (-0.5*npar*log(2*pi)), rep(0,npar), rep(-1,npar))

    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
        length(aa.start) ), ndeps = rep.int(0.001, length(aa.start)), maxit = 100L, abstol = -Inf, 
        reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
        gamma = 2, REPORT = 10, type = 1, lmm = 5, factr = 1e+07, 
        pgtol = 0, tmax = 10, temp = 10)  
        
    con[(namc <- names(control2))] <- control2            
    result <- mcll_coeff(start=aa.start, data.t=data.t,  par.t= vec.t, alp=alp,  
                    nodes=nodes, weights=weights, method=method2, control=con )  

    # return values of mcll_coeff 
    pol_coeff <- as.vector(result$a_coeff)
    convergence <- result$convergence
    value <- result$value
    counts <- result$counts
    message <- result$message 
                    
    linear <-sum(t(vec.t)*pol_coeff[2:(npar+1)])
    quad <- sum(t(vec.t)^2*pol_coeff[(1+npar+1): length(pol_coeff)]*0.5 )
    pol <- pol_coeff[1] + linear + quad 

    log.post <- pol 
        
    # prior 
    log.prior <- prior.func(vec)
        
    # log likelihood 
    log.lik <- log.post - log.prior -log(det(ch))

    return(-log.lik)  
}
