bbeta <- function(formula, 
                  burn.in = 100, 
                  main.iters = 1000, 
                  m.prior = NULL, 
                  sigma.prior = NULL, 
                  sigma.epsilon = NULL,
                  ...){
                   	
    y <- ergm.getnetwork(formula)
    model <- ergm.getmodel(formula,y)
    stats0 <- summary(formula)
    nstats <- length(stats0)
    
    if (is.null(m.prior)) 
        m.prior <- rep(0, nstats)
    if (is.null(sigma.prior)) 
        sigma.prior <- diag(100, nstats)
	if (is.null(sigma.epsilon)) 
    	sigma.epsilon <- diag(0.025, nstats)
    
    Theta <- matrix(NA, nrow = main.iters, ncol = nstats)
    theta <- runif(nstats, min = -0.1, max = 0.1)
    acc.counts <- 0 
    tot.iters <- burn.in + main.iters

    for (k in 1:tot.iters) {

            theta1 <- rmvnorm(1, mean = theta, sigma = sigma.epsilon)[1, ]
            
            pr <- dmvnorm(rbind(theta1, theta), mean = m.prior, sigma = sigma.prior, log = TRUE)

            beta <- (theta1 - theta) %*% stats0 - 
                    sum(log(1 + exp(apply(combn(theta1, 2), 2, sum)))) +
                    sum(log(1 + exp(apply(combn(theta, 2), 2, sum)))) +
                    pr[1] - pr[2]
         
            if (beta > log(runif(1))) {
                theta <- theta1
                if (k > burn.in) 
                  acc.counts <- acc.counts + 1
            }
    
    Theta[k - burn.in,] <- theta        
	}

    out = list(
        formula = formula,
        model = model,
        specs = model$coef.names,
        stats = stats0, 
        Theta = Theta,
        acc.rate = acc.counts/main.iters,
        nstats = nstats)

    out
}
