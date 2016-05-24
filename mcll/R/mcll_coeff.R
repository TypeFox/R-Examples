# Standard error estimation
mcll_coeff <- function(start, data.t, par.t,  alp,  nodes=nodes, weights=weights, 
            method,  lower, upper, control ) { 

    # data.t: transformed posterior 
    # par.t: transformed parameter values 
    # alp: alpha smoothing parameter 
    # nodes: quadrature locations 
    # weights: quadrature weights at the nodes 
    # method: optimization method 
    # lower, upper: boundary values  for the "L-BFGS-B" method
    # control: control options 

    data.t <- as.matrix(data.t)
    
    # optim options 
       
    if (method=="L-BFGS-B") {
        result <-optim(par=start, fn=lik_local, x=par.t, data.t=data.t, alp= alp, 
            nodes=nodes, weights=weights, low=-3, up=3, 
            method=method,  lower = lower, upper = upper, control=control  )
    } else {
        result <-optim(par=start, fn=lik_local, x=par.t, data.t=data.t, alp= alp, 
            nodes=nodes, weights=weights, low=-3, up=3, 
            method=method,  control=control  )
    }
    

    # return 
    a_coeff <- result$par
    convergence <- result$convergence
    value <- result$value
    counts <- result$counts
    message <- result$message 
    
    
    return(list("a_coeff"=a_coeff, "convergence"=convergence, "value"=value, "counts"=counts, "message"=message)) 

}
