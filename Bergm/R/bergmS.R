bergmS <- function(formulae,
                   iters=10000,
                   m.priors=NULL,
                   sigma.priors=NULL,
                   gammas=NULL,
                   nchains=NULL,
                   sigma.epsilons=NULL,
                   aux.iters=1000,
                   main.iters=NULL,
                   burn.ins=NULL,
                   ...){
   
    ptm = proc.time()
    nmodels <- length(formulae)
    control <- control.simulate.formula(MCMC.burnin = aux.iters, 
                                        MCMC.interval = 0)
    control$MCMC.samplesize <- 1
    dims <- rep(0, nmodels)
    y <- ergm.getnetwork(formulae[[1]])
    
    if (is.null(main.iters)) 
        main.iters <- rep(1000, nmodels)
    if (is.null(burn.ins)) 
        burn.ins <- rep(100, nmodels)
    if (is.null(gammas)) 
        gammas <- rep(0.5, nmodels)
    
    models <- vector("list", nmodels)
    specs <- vector("list", nmodels) 
    theta <- vector("list", nmodels) 
    Sigma <- vector("list", nmodels) 
    Mu <- vector("list", nmodels) 
    post <- vector("list", nmodels) 
    stats <- vector("list", nmodels)

    for (i in 1L:nmodels) {
        post[[i]] <- bergm(formula = formulae[[i]], burn.in = burn.ins[i], 
            main.iters = main.iters[i], aux.iters = aux.iters, 
            m.prior = m.priors[[i]], sigma.prior = sigma.priors[[i]], 
            gamma = gammas[i], nchains = nchains[i], sigma.epsilon = sigma.epsilons[[i]])
        if (post[[i]]$nchains > 1) {
            Sigma[[i]] <- cov(apply(post[[i]]$Theta, 2, cbind))
            diag(Sigma[[i]]) = 8/7 * diag(Sigma[[i]])
        }
        else {
            Sigma[[i]] <- var(post[[i]]$Theta)
            diag(Sigma[[i]]) = Sigma[[i]]
        }
        Mu[[i]] <- apply(post[[i]]$Theta, 2, mean)
        theta[[i]] <- runif(post[[i]]$dim, min = -0.1, max = 0.1)
        dims[i] <- post[[i]]$dim
        specs[[i]] <- post[[i]]$specs
        
        #bergm.output(post[[i]])
    }
    
    #
    
        Theta <- matrix(0L, iters, max(dims))
        M <- rep.int(0, iters)
        m <- sample(seq(1:nmodels), 1)
        Wrate <- matrix(0L, 2, nmodels)
        Brate <- matrix(0L, 2, 1)

        for (i in 1L:iters) {
            m1 <- sample(1:nmodels, 1)
            if (m1 != m) {
                Brate[2, ] <- Brate[2, ] + 1
            }
            else {
                Wrate[2, m] <- Wrate[2, m] + 1
            }
            theta1 <- rmvnorm(1, mean = Mu[[m1]], sigma = Sigma[[m1]])[1,]
            ww1 <- dmvnorm(theta1, mean = Mu[[m1]], 
                           sigma = Sigma[[m1]],log=TRUE)
            pr1 <- dmvnorm(theta1, mean = post[[m1]]$m.prior, 
                           sigma = post[[m1]]$sigma.prior,log=TRUE)
            ww <- dmvnorm(theta[[m]], mean = Mu[[m]], 
                          sigma = Sigma[[m]],log=TRUE)
            pr <- dmvnorm(theta[[m]], mean = post[[m]]$m.prior, 
                          sigma = post[[m]]$sigma.prior,log=TRUE)
                          
            z <- ergm.mcmcslave(Clist = post[[m1]]$Clist, MHproposal = post[[m1]]$MHproposal, 
                                eta0 = theta1, control = post[[m1]]$control, verbose = FALSE)
            
            y1 <- newnw.extract(y, z)
            smy1 <- ergm.getglobalstats(y1, post[[m]]$model)
           
            delta <- smy1 - post[[m]]$stats
            beta <- (theta[[m]] %*% delta) + 
                    (theta1 %*% as.numeric(-z$s)) + 
                    pr1 - pr + ww - ww1
                    
            if (beta >= log(runif(1))) {
                theta[[m1]] <- theta1
                if (m1 != m) {
                    Brate[1, ] <- Brate[1, ] + 1
                    m <- m1
                }
                else {
                    Wrate[1, m] <- Wrate[1, m] + 1
                }
            }
            Theta[i, 1:dims[m]] <- theta[[m]]
            M[i] <- m
        }
        time = proc.time() - ptm
        out = list(M = M, iters = iters, Theta = Theta, post = post, formulae = formulae, 
            models = models, specs = specs, dims = dims, nmodels = nmodels, 
            Baccept = Brate[1, ]/Brate[2, ], Waccept = Wrate[1, ]/Wrate[2, 
                ], time = time)
        out

}