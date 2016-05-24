bergm <- function (formula, 
                   burn.in=10,
                   main.iters=1000,
                   aux.iters=1000, 
                   m.prior = NULL, 
                   sigma.prior = NULL, 
                   nchains = NULL, 
                   gamma = 0.5, 
                   sigma.epsilon = NULL,
                   ...){ 	
                   	
    y <- ergm.getnetwork(formula)
    model <- ergm.getmodel(formula, y)
    Clist <- ergm.Cprepare(y, model)
        
    stats0 <- summary(formula)
    control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                        MCMC.interval = 0)
    control$MCMC.samplesize <- 1
    
    MHproposal <- MHproposal.ergm(object= model, 
                                  constraints = ~., arguments = control$MCMC.prop.args, 
                                  nw = y, weights = control$MCMC.prop.weights, class = "c", 
                                  reference = ~Bernoulli, response = NULL)     
    
    snooker <- 0
    if (is.null(m.prior)) m.prior <- rep(0, Clist$nstats)
    if (is.null(sigma.prior)) sigma.prior <- diag(100, Clist$nstats)
    if (is.null(nchains)) nchains <- 2 * Clist$nstats
    if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, Clist$nstats)
    if (Clist$nstats == 1) {
    	 nchains <- 1
        sigma.epsilon <- diag(gamma, Clist$nstats)
    }
    Theta <- array(NA, c(main.iters, Clist$nstats, nchains))
    # MPLE as starting point
    theta <- matrix(runif(Clist$nstats * nchains, min = -0.001, max = 0.001), Clist$nstats, nchains) +
             matrix(rep(ergm(formula, estimate = "MPLE")[[1]],nchains),Clist$nstats,nchains)
    acc.counts <- rep(0L, nchains)
    theta1 <- rep(NA, Clist$nstats)
    tot.iters <- burn.in + main.iters

    for (k in 1L:tot.iters) {
        for (h in 1L:nchains) {
            if (Clist$nstats > 1 && nchains > 1) {
                snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
            }
            
            theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1,]
            
            pr <- dmvnorm(rbind(theta1, theta[, h]), 
                          mean = m.prior, 
                          sigma = sigma.prior,
                          log=TRUE)
                          
            delta <- as.numeric(ergm.mcmcslave(Clist, 
                                    MHproposal, 
                                    eta0 = theta1, 
                                    control, 
                                    verbose = FALSE)$s)
                                    
            beta <- (theta[, h] - theta1) %*% delta + pr[1] - pr[2]
            
            if (beta >= log(runif(1))) {
                
                theta[, h] <- theta1
                if (k > burn.in) acc.counts[h] <- acc.counts[h] + 1
            
            }
        
        }
        if (k > burn.in) Theta[k - burn.in, , ] <- theta
    }
    if (nchains == 1) Theta <- as.matrix(Theta[, , 1])


    out = list(Clist = Clist, MHproposal = MHproposal, control = control, 
        formula = formula, model = model, nnodes = Clist$n, specs = model$coef.names, 
        dim = Clist$nstats, nchains = nchains, stats = stats0, 
        Theta = Theta, nchains = nchains, acc.rates = acc.counts / main.iters, 
        m.prior = m.prior, sigma.prior = sigma.prior, aux.iters = aux.iters)
    out

}
