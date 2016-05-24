topt <- function(
    x, #initial support points 
    w, #initial weights
    eta, #list of regression models
    theta.fix, #list of fixed parameters
    theta.var, #initial values of variable parameters
    p, #comparison table
    x.lb, #left border for support points
    x.rb, #right border for support points
    opt #other options
    )
{
    ### Initialization ###
    
    #Options
    method <- ifelse(is.null(opt$method), 1, opt$method)
    max.iter <- ifelse(is.null(opt$max.iter), 100, opt$max.iter)
    des.eff <- ifelse(is.null(opt$des.eff), 0.999, opt$des.eff)
    derivative.epsilon <- ifelse(is.null(opt$derivative.epsilon), 1e-5, opt$derivative.epsilon)
    support.epsilon <- ifelse(is.null(opt$support.epsilon), .Machine$double.eps^0.25, opt$support.epsilon)
    weights.evaluation.epsilon <- ifelse(is.null(opt$weights.evaluation.epsilon), ifelse(method == 1, 1e-3, 1e-7), opt$weights.evaluation.epsilon)
    weights.evaluation.max.iter <- ifelse(is.null(opt$weights.evaluation.max.iter), ifelse(method == 1, 10, 100), opt$weights.evaluation.max.iter)

    start.time <- proc.time()

    efficiency <- numeric(max.iter)
    functional <- numeric(max.iter)
    theta.var <- array(list(NULL), c(length(eta),length(eta)))
    if(method == 2)
    {
        epsilon <- array(list(NULL), c(length(eta),length(eta)))
        jacob <- array(list(NULL), c(length(eta),length(eta)))
    }
    for(i in 1:length(eta))
        for(j in 1:length(eta))
            if(p[i,j] != 0)
            {
                theta.var[[i,j]] <- optim(
                    par = theta.fix[[j]],
                    function(theta) Tfs(x, w, eta[[i]], eta[[j]], theta.fix[[i]], theta)
                    )$par
            }

    ### Main loop ###

    iter <- 0 
    delta.0 <- 0
    psi.0 <- 1
    while(iter < max.iter)
    {
        iter <- iter + 1

        local.maximums <- evaluate.local.maximums(x, derivative.epsilon, p, eta, theta.fix, theta.var, x.lb, x.rb)

        #Efficiency estimation

        s <- local.maximums[which.max(Psi(local.maximums, p, eta, theta.fix, theta.var))]
        delta.0 <- Tp(x, w, p, eta, theta.fix, theta.var)
        psi.0 <- Psi(s, p, eta, theta.fix, theta.var)
        efficiency[iter] <- delta.0 / psi.0
        functional[iter] <- delta.0
        if(delta.0 / psi.0 >= des.eff)
            break

        x <- unique(sort(c(x, local.maximums)))

        tmp <- switch(
            method,
            evaluate.weights.1(x, p, eta, theta.fix, theta.var, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon),
            evaluate.weights.2(x, p, eta, theta.fix, theta.var, epsilon, jacob, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon)
        )
        x <- tmp$x
        w <- tmp$w
        theta.var <- tmp$theta.var
    }
    x <- sort(local.maximums)
    tmp <- switch(
        method,
        evaluate.weights.1(x, p, eta, theta.fix, theta.var, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon),
        evaluate.weights.2(x, p, eta, theta.fix, theta.var, epsilon, jacob, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon)
    )
    x <- tmp$x
    w <- tmp$w
    theta.var <- tmp$theta.var
    
    end.time <- proc.time()
    
    list(
        x = x, 
        w = w, 
        efficiency = efficiency[efficiency > 0],
        functional = functional[functional > 0],
        eta = eta,
        theta.fix = theta.fix,
        theta.var = theta.var,
        p = p,
        x.lb = x.lb,
        x.rb = x.rb,
        max.iter = max.iter,
        done.iter = iter,
        des.eff = des.eff,
        time = end.time - start.time
    )
}
