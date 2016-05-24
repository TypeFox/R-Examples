# Parameter estimation 

mcll_est <- function(data, prior.func, alp=0.7, 
        method="BFGS",  lower = -Inf, upper = Inf, 
        control=list() , use.locfit=TRUE, con.manual=list(method="BFGS", lower = -Inf, upper = Inf, control=list()) ) { 

    # data: posterior samples 
    # prior.func: prior function (user specified) that provies a log prior density 
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

    # prior function
    prior.func <- match.fun(prior.func)
    
    # posterior mean as starting values 
    par.start <- apply(as.matrix(data.t), 2, mean) 
 
    # optim option
    con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
        length(bb) ), ndeps = rep.int(0.001, length(bb)), maxit = 100L, abstol = -Inf, 
        reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
        gamma = 2, REPORT = 10, type = 1, lmm = 5, factr = 1e+07, 
        pgtol = 0, tmax = 10, temp = 10)

    con[(namc <- names(control))] <- control 
        
    if (use.locfit==TRUE) {

        if (method=="L-BFGS-B") {
            try(result <- optim(par=par.start, fn=lik_locfit, data.t=data.t, ch=ch, bb=bb,  prior.func= prior.func, alp=alp, 
                        method=method,  lower = lower, upper = upper, control=con)
                         , silent=TRUE)
        } else {            
            try(result <- optim(par=par.start, fn=lik_locfit, data.t=data.t, ch=ch, bb=bb, prior.func= prior.func, alp=alp, 
                        method=method,  control=con)
                         , silent=TRUE)
        }  
        
    } else {
    
        # starting values 
        npar <- length(par.start)  # dimension 
        aa <- c( (-0.5*npar*log(2*pi)), rep(0,npar), rep(-1,npar))
    
        aa.con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
            length(aa) ), ndeps = rep.int(0.001, length(aa)), maxit = 100L, abstol = -Inf, 
            reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
            gamma = 2, REPORT = 10, type = 1, lmm = 5, factr = 1e+07, 
            pgtol = 0, tmax = 10, temp = 10)  
        
        aa.con[(namc <- names(con.manual$control))] <- con.manual$control
      
        if (method=="L-BFGS-B") {
            try(result <- optim(par=par.start, fn=lik_manual, data.t=data.t, ch=ch, bb=bb,  prior.func= prior.func, alp=alp, 
                        method=method,  lower = lower, upper = upper, control=con, 
                        method2=con.manual$method,  lower2 = con.manual$lower, upper2 = con.manual$upper, control2=aa.con)
                         , silent=TRUE)
        } else {      
            try(result <- optim(par=par.start, fn=lik_manual, data.t=data.t, ch=ch, bb=bb,  prior.func= prior.func, alp=alp, 
                        method=method,  control=con, 
                        method2=con.manual$method,  lower2 = con.manual$lower, upper2 = con.manual$upper, control2=aa.con)
                         , silent=TRUE)
        }    
  
    }
    
    if (is.na(result[2])==FALSE) {   
        # parameter estimates               
        par0 <- result$par
        # transformation back      
        par.t <- t((ch%*%par0) +bb)  
        
        # optimization return values 
        convergence <- result$convergence  
        value <- result$value
        counts <- result$counts
        message <- result$message   
      
      
        return(list("par" = par.t, "convergence"=convergence, "value"=value, "counts"=counts, "message"=message)) 
    } else {
    
        print("Error during optimization")
    }
    
    
}
