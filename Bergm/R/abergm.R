abergm <- function (formula, 
                    burn.in = 100, 
                    main.iters = 1000, 
                    aux.iters = 1000, 
                    m.prior = NULL, 
                    sigma.prior = NULL, 
                    nchains = NULL, 
                    gamma = 0.5, 
                    method = 'Adaptive.chains', 
                    rectangular = TRUE,
                    sigma.epsilon = NULL, 
                    updategap = 10,
                    ...){
    	
	methodnames <- c('ADS','Adaptive.past','Adaptive.chains')
	stopifnot(method %in% methodnames)
	
	UpdateCov<- function(tv,oldcov, d, t){
		sdt <- 2.4^2/(d*t)
		xb <- colSums(tv)
		x <- tv[t,]
		xob <- xb - x
		epsilon <- .1^20
		newcov <- (t - 1)/t * oldcov + sdt * (tcrossprod(xob)/t - tcrossprod(xb)/(t+1) + tcrossprod(x) + epsilon * diag(d))
	}	
	
    y <- ergm.getnetwork(formula)
    model <- ergm.getmodel(formula, y)
    Clist <- ergm.Cprepare(y, model)
    stats0 <- summary(formula)
    control <- control.simulate.formula(MCMC.burnin = aux.iters, MCMC.interval = 0)
    control$MCMC.samplesize <- 1
    MHproposal <- MHproposal.ergm(object = model, constraints = ~., 
        arguments = control$MCMC.prop.args, nw = y, weights = control$MCMC.prop.weights, 
        class = "c", reference = ~Bernoulli, response = NULL)

    if (is.null(m.prior)) 
        m.prior <- rep(0, Clist$nstats)
    if (is.null(sigma.prior)) 
        sigma.prior <- diag(100, Clist$nstats)
    if (is.null(nchains)) 
        nchains <- 2 * Clist$nstats
    if (is.null(sigma.epsilon)) 
        sigma.epsilon <- diag(0.0025, Clist$nstats)
    if (Clist$nstats == 1) {
        nchains <- 1
        sigma.epsilon <- diag(gamma, Clist$nstats)
    }
    Theta <- array(NA, c(main.iters, Clist$nstats, nchains))
    theta <- matrix(runif(Clist$nstats * nchains, min = -0.1, max = 0.1), Clist$nstats, nchains)
    acc.counts <- rep(0, nchains)
    acc.counts2 <- rep(0, nchains)
    tot.iters <- burn.in + main.iters
    iters2 <- 0 
    
    SEQ <- (burn.in + 50):tot.iters
    COV <- NA
    snooker <- 0
    updatesteps <- seq(1,tot.iters,updategap)

    for (k in 1L:tot.iters) {
        for (h in 1L:nchains) {
            if (Clist$nstats > 1 & nchains > 1) {

            	if(k >= (burn.in+50) & method %in% c('Adaptive.past','Adaptive.chains') & k %in% updatesteps){   
            		if(method=='Adaptive.past'){ 
            			if(any(is.na(COV))){ 
            				if(rectangular)
            					COV <- cov(apply(Theta[1:(k-burn.in-1),,],2,rbind))
            				else             
            					COV <- cov(Theta[1:(k-burn.in-1),,h])              
            			}else{ 
            				if(rectangular)
            					tv <- apply(Theta[1:(k-burn.in-1),,],2,rbind) 
            				else  
            					tv <- Theta[1:(k-burn.in-1),,h]                
            				
            				COV <- UpdateCov(tv,COV,Clist$nstats,nrow(tv))
            			}
            			SIGMA <- COV*((2.38)^2/Clist$nstats)

            			FT <- sample(c(TRUE,FALSE),1,prob=c(.99,.01))
            			sigma.epsilon <- FT*SIGMA + (1-FT)*diag(0.0025, Clist$nstats)
         			
            		}else{
            			if(method=='Adaptive.chains')   		
            				sigma.epsilon <- (1/Clist$nstats)*cov(t(theta))*((2.38)^2)
            		}
            	}
            	if(k > burn.in | method %in% c('Adaptive.past','Adaptive.chains')){
            		snooker <- 0
            	}
            	if(method=='ADS')
                		snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff) 
            }
            err <- rmvnorm(1, sigma = sigma.epsilon)[1,] 
            theta1 <- theta[,h] + snooker + err
            delta1 <- ergm.mcmcslave(Clist, MHproposal, eta0 = theta1, control, verbose = FALSE)$s    
            pr1 <- diff(dmvnorm(rbind(theta[,h],theta1),mean = m.prior,sigma=sigma.prior,log=TRUE))
            beta1 <- (theta[,h]-theta1) %*% delta1 + pr1
            
            if(beta1 >= log(runif(1))){
                theta[,h] <- theta1
                if (k > burn.in)
                		acc.counts[h] <- acc.counts[h] + 1
            }
          	if(beta1 < log(runif(1)) & k > burn.in){
            	iters2 <- iters2+1
            	theta2 <- theta[,h] - snooker + rmvnorm(1, sigma = sigma.epsilon/2)[1,]
            	pr1.num <- diff(dmvnorm(rbind(theta2,theta1),mean = m.prior,sigma = sigma.prior,log=TRUE))
            	beta1.num <- (theta2-theta1) %*% delta1 + pr1.num
				if(beta1.num < 0){
					t2t0 <- rbind(theta[,h],theta2)
            		beta2 <- (theta[,h]-theta2) %*% 
            		         ergm.mcmcslave(Clist, MHproposal, eta0 = theta2, control, verbose = FALSE)$s + 
            		         pr1 - pr1.num + 
            		         diff(dmvnorm(t2t0, mean = theta1-snooker, sigma = sigma.epsilon,log=TRUE)) + 
            		         log((1-exp(beta1.num))/(1-exp(beta1)))
          		 
            			if(beta2 >= log(runif(1))){
                        theta[,h] <- theta2
                    		if(k > burn.in) 
                    			acc.counts2[h] <- acc.counts2[h] + 1
                	}
            		}
          	}
        }
        if (k > burn.in) 
            Theta[k - burn.in, , ] <- theta
    }
    if (nchains == 1)
        Theta <- as.matrix(Theta[, , 1])
    
    out = list(Clist = Clist, MHproposal = MHproposal, control = control, 
        formula = formula, model = model, nnodes = Clist$n, specs = model$coef.names, 
        dim = Clist$nstats, nchains = nchains, stats = stats0, 
        Theta = Theta, nchains = nchains, acc.rate = acc.counts/main.iters, 
        m.prior = m.prior, sigma.prior = sigma.prior, aux.iters = aux.iters, 
        acc.rate2 = acc.counts2/iters2)
    out
}
