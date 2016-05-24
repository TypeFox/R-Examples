binomial.CARlocalised <- function(formula, data=NULL, G, trials,  W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.delta=NULL, prior.tau2=NULL, verbose=TRUE)
    {
    #### Check on the verbose option
    if(is.null(verbose)) verbose=TRUE     
    if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)
    
    if(verbose)
    {
        cat("Setting up the model\n")
        a<-proc.time()
    }else{}
    
    
    
    ##############################################
    #### Format the arguments and check for errors
    ##############################################
    #### Overall formula object
    frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
    if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths or the data object has not been specified.", call.=FALSE)
       
    #### Format and check the neighbourhood matrix W
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    K <- nrow(W)
    if(ncol(W)!= K) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
    
    #### Response variable
    Y <- model.response(frame)
    N.all <- length(Y)
    N <- N.all / K
    if(floor(N.all/K)!=ceiling(N.all/K)) stop("The number of spatial areas is not a multiple of the number of data points.", call.=FALSE)
    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
    int.check <- N.all-sum(ceiling(Y)==floor(Y))
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)   
    failures <- trials - Y
    Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE) 
    failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=FALSE)  
    which.miss <- as.numeric(!is.na(Y))
    which.miss.mat <- matrix(which.miss, nrow=K, ncol=N, byrow=FALSE)
     
    #### Offset variable
    ## Create the offset
    offset <- try(model.offset(frame), silent=TRUE)
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,N.all)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
    
    
    #### Design matrix
    ## Create the matrix
    X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
    if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
    if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
    ptemp <- ncol(X)
    
    if(ptemp==1)
    {
        X <- NULL
        regression.vec <- rep(0, N.all)
        regression.mat <- matrix(regression.vec, nrow=K, ncol=N, byrow=FALSE)
        p <- 0
    }else
    {
        ## Check for linearly related columns
        cor.X <- suppressWarnings(cor(X))
        diag(cor.X) <- 0
        if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
        if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
        if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
        
        ## Remove the intercept term
        int.which <- which(apply(X,2,sd)==0)
        colnames.X <- colnames(X)
        X <- as.matrix(X[ ,-int.which])
        colnames(X) <- colnames.X[-int.which]
        p <- ncol(X)
        
        ## Standardise X
        X.standardised <- X
        X.sd <- apply(X, 2, sd)
        X.mean <- apply(X, 2, mean)
        X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
        
        for(j in 1:p)
        {
            if(length(table(X[ ,j]))>2)
            {
                X.indicator[j] <- 1
                X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
            }else
            {
                X.indicator[j] <- 0
            }
        }
        
        ## Compute a starting value for beta
        dat <- cbind(Y, failures)
        mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
        beta.mean <- mod.glm$coefficients
        beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
        beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
        regression.vec <- X.standardised %*% beta
        regression.mat <- matrix(regression.vec, nrow=K, ncol=N, byrow=FALSE)   
    }
    
    
    #### Format and check the number of clusters G     
    if(length(G)!=1) stop("G is the wrong length.", call.=FALSE)    
    if(!is.numeric(G)) stop("G is not numeric.", call.=FALSE)    
    if(G<=1) stop("G is less than 2.", call.=FALSE)    
    if(G!=round(G)) stop("G is not an integer.", call.=FALSE) 
    if(floor(G/2)==ceiling(G/2))
    {
        Gstar <- G/2    
    }else
    {
        Gstar <- (G+1)/2          
    }
    
    
    #### Format and check the MCMC quantities
    if(is.null(burnin)) stop("the burnin argument is missing", call.=FALSE)
    if(is.null(n.sample)) stop("the n.sample argument is missing", call.=FALSE)
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
    if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
    if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 
    
    
    #### Check and specify the priors
    if(!is.null(X))
    {
        if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
        if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
        if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
        if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)       
        
        if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
        if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
        if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
        if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
        if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE) 
    }else
    {}
    
    if(is.null(prior.delta)) prior.delta <- 10       
    if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
    if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    
    
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
    if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
    
    
    
    #### Specify the initial parameter values
    theta.hat <- Y / trials
    theta.hat[theta.hat==0] <- 0.01
    theta.hat[theta.hat==1] <- 0.99
    res.temp <- log(theta.hat / (1 - theta.hat)) - regression.vec - offset
    res.sd <- sd(res.temp, na.rm=TRUE)/5
    phi.mat <- matrix(rnorm(n=N.all, mean=0, sd = res.sd), nrow=K, byrow=FALSE)
    tau2 <- var(phi)/10
    gamma <- runif(1)
    Z <- sample(1:G, size=N.all, replace=TRUE)
    Z.mat <- matrix(Z, nrow=K, ncol=N, byrow=FALSE)
    lambda <- sort(runif(G, min=min(res.temp), max=max(res.temp)))
    lambda.mat <- matrix(rep(lambda, N), nrow=N, byrow=TRUE)
    delta <- runif(1,1, min(2, prior.delta))
    mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)


    
    ## Compute the blocking structure for beta     
    if(!is.null(X))
    {
        blocksize.beta <- 5 
        if(blocksize.beta >= p)
        {
            n.beta.block <- 1
            beta.beg <- 1
            beta.fin <- p
        }else
        {
            n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
            remainder <- p - n.standard * blocksize.beta
            
            if(remainder==0)
            {
                beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
                beta.fin <- seq(blocksize.beta, p, blocksize.beta)
                n.beta.block <- length(beta.beg)
            }else
            {
                beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
                beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
                n.beta.block <- length(beta.beg)
            }
        }
    }else{}
    
    
    #### Set up matrices to store samples
    n.keep <- floor((n.sample - burnin)/thin)
    samples.Z <- array(NA, c(n.keep, N.all))
    samples.lambda <- array(NA, c(n.keep, G))
    samples.delta <- array(NA, c(n.keep, 1))
    samples.tau2 <- array(NA, c(n.keep, 1))
    samples.gamma <- array(NA, c(n.keep, 1))
    samples.phi <- array(NA, c(n.keep, N.all))
    samples.fitted <- array(NA, c(n.keep, N.all))
    samples.deviance <- array(NA, c(n.keep, 1))
    samples.like <- array(NA, c(n.keep, N.all))
    
    if(!is.null(X))
    {
        samples.beta <- array(NA, c(n.keep, p))
        accept.all <- rep(0,8)
        proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
        chol.proposal.corr.beta <- chol(proposal.corr.beta) 
        proposal.sd.beta <- 0.01
    }else
    {
        accept.all <- rep(0,6)    
    }
    
    accept <- accept.all
    proposal.sd.lambda <- 0.1
    proposal.sd.delta <- 0.1
    proposal.sd.phi <- 0.1
    Y.extend <- matrix(rep(Y, G), byrow=F, ncol=G)
    delta.update <- matrix(rep(1:G, N.all-K), ncol=G, byrow=T)
    tau2.posterior.shape <- prior.tau2[1] + N * (K-1) /2
    
     
    
    #### Spatial quantities
    ## Create the triplet object
    W.triplet <- c(NA, NA, NA)
    for(i in 1:K)
    {
        for(j in 1:K)
        {
            if(W[i,j]>0)
            {
                W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
            }else{}
        }
    }
    W.triplet <- W.triplet[-1, ]     
    W.n.triplet <- nrow(W.triplet) 
    W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
    W.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
    
    
    ## Create the start and finish points for W updating
    W.begfin <- array(NA, c(K, 2))     
    temp <- 1
    for(i in 1:K)
    {
        W.begfin[i, ] <- c(temp, (temp + W.neighbours[i]-1))
        temp <- temp + W.neighbours[i]
    }
    
    
    
    
    ###########################
    #### Run the Bayesian model
    ###########################
    ## Start timer
    if(verbose)
    {
        cat("Generating", n.sample, "samples\n", sep = " ")
        progressBar <- txtProgressBar(style = 3)
        percentage.points<-round((1:100/100)*n.sample)
    }else
    {
        percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
    
    for(j in 1:n.sample)
    {
        ####################
        ## Sample from beta
        ####################
        if(!is.null(X))
        {
            proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
            proposal.beta <- beta
            offset.temp <- offset + as.numeric(mu) + as.numeric(phi.mat)   
            for(r in 1:n.beta.block)
            {
                proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
                prob <- binomialbetaupdate(X.standardised, N.all, p, beta, proposal.beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta, which.miss)
                if(prob > runif(1))
                {
                    beta[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
                    accept[7] <- accept[7] + 1  
                }else
                {
                    proposal.beta[beta.beg[r]:beta.fin[r]] <- beta[beta.beg[r]:beta.fin[r]]
                }
            }
            accept[8] <- accept[8] + n.beta.block    
            regression.vec <- X.standardised %*% beta
            regression.mat <- matrix(regression.vec, nrow=K, ncol=N, byrow=FALSE)   
        }else{}
        
        
        
        #######################     
        #### Sample from lambda
        #######################
        #### Propose a new value
        proposal.extend <- c(-100, lambda, 100) 
        for(r in 1:G)
        {
            proposal.extend[(r+1)] <- rtrunc(n=1, spec="norm", a=proposal.extend[r], b=proposal.extend[(r+2)], mean=proposal.extend[(r+1)], sd=proposal.sd.lambda)
        }
        proposal <- proposal.extend[-c(1, (G+2))]
        
        #### Compute the data likelihood
        lp.current <- lambda[Z] + offset + as.numeric(regression.mat) + as.numeric(phi.mat)
        lp.proposal <- proposal[Z] + offset + as.numeric(regression.mat) + as.numeric(phi.mat)
        p.current <- exp(lp.current) / (1 + exp(lp.current))
        p.proposal <- exp(lp.proposal) / (1 + exp(lp.proposal))
        like.current <- Y * log(p.current) + failures * log(1 - p.current)
        like.proposal <- Y * log(p.proposal) + failures * log(1 - p.proposal)
        prob <- exp(sum(like.proposal - like.current))
        if(prob > runif(1))
        {
            lambda <- proposal
            lambda.mat <- matrix(rep(lambda, N), nrow=N, byrow=TRUE)
            mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
            accept[1] <- accept[1] + 1  
        }else
        {
        }
        accept[2] <- accept[2] + 1           
        
        
        
        ##################     
        #### Sample from Z
        ##################
        prior.offset <- rep(NA, G)
        for(r in 1:G)
        {
            prior.offset[r] <-  log(sum(exp(-delta * ((1:G - r)^2 + (1:G - Gstar)^2)))) 
        }
        mu.offset <- offset.mat + regression.mat + phi.mat
        test <- Zupdatesqbin(Z=Z.mat, Offset=mu.offset, Y=Y.mat, delta=delta, lambda=lambda, nsites=K, ntime=N, G=G, SS=1:G, prioroffset=prior.offset, Gstar=Gstar, failures=failures.mat)          
        Z.mat <- test
        Z <- as.numeric(Z.mat)
        mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
        
        
        
        ######################
        #### Sample from delta
        ######################
        proposal.delta <-  rtrunc(n=1, spec="norm", a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)
        sum.delta1 <- sum((Z - Gstar)^2)
        sum.delta2 <- sum((Z.mat[ ,-1] - Z.mat[ ,-N])^2)
        current.fc1 <- -delta * (sum.delta1 + sum.delta2) - K *  log(sum(exp(-delta * (1:G - Gstar)^2))) 
        proposal.fc1 <- -proposal.delta * (sum.delta1 + sum.delta2) - K *  log(sum(exp(-proposal.delta * (1:G - Gstar)^2)))                 
        Z.temp <- matrix(rep(as.numeric(Z.mat[ ,-N]),G), ncol=G, byrow=FALSE)
        Z.temp2 <- (delta.update - Z.temp)^2 + (delta.update - Gstar)^2
        current.fc <- current.fc1 - sum(log(apply(exp(-delta * Z.temp2),1,sum)))
        proposal.fc <- proposal.fc1 - sum(log(apply(exp(-proposal.delta * Z.temp2),1,sum)))        
        prob <- exp(proposal.fc - current.fc)       
        if(prob > runif(1))
        {
            delta <- proposal.delta
            accept[3] <- accept[3] + 1  
        }else
        {
        }
        accept[4] <- accept[4] + 1   
        
        
        
        ####################
        #### Sample from phi
        ####################
        phi.offset <- mu + offset.mat + regression.mat
        temp1 <- binomialarcarupdate(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, gamma, 1, Y.mat, failures.mat, proposal.sd.phi, phi.offset, W.triplet.sum, which.miss.mat)      
        phi.temp <- temp1[[1]]
        phi <- as.numeric(phi.temp)
            for(i in 1:G)
            {
            phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
            }
        
        phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
        accept[5] <- accept[5] + temp1[[2]]
        accept[6] <- accept[6] + K*N    
  
        
        
        ####################
        ## Sample from gamma
        ####################
        temp2 <- gammaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, 1)
        mean.gamma <- temp2[[1]] / temp2[[2]]
        sd.gamma <- sqrt(tau2 / temp2[[2]]) 
        gamma <- rtrunc(n=1, spec="norm", a=0, b=1, mean=mean.gamma, sd=sd.gamma)   

        
        
        ####################
        ## Samples from tau2
        ####################
        temp3 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, 1, gamma)
        tau2.posterior.scale <- temp3 + prior.tau2[2] 
        tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))          
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        lp <- as.numeric(mu + offset.mat + regression.mat + phi.mat)
        prob <- exp(lp) / (1+exp(lp))
        fitted <- trials * prob
        deviance.all <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)
        like <- exp(deviance.all)
        deviance <- -2 * sum(deviance.all) 
        
        
        ###################
        ## Save the results
        ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
            ele <- (j - burnin) / thin
            samples.delta[ele, ] <- delta
            samples.lambda[ele, ] <- lambda
            samples.Z[ele, ] <- Z
            samples.phi[ele, ] <- as.numeric(phi.mat)
            samples.tau2[ele, ] <- tau2
            samples.gamma[ele, ] <- gamma
            samples.deviance[ele, ] <- deviance
            samples.fitted[ele, ] <- fitted
            samples.like[ele, ] <- like
            if(!is.null(X)) samples.beta[ele, ] <- beta        
        }else
        {
        }
        
        
        ########################################
        ## Self tune the acceptance probabilties
        ########################################
        k <- j/100
        if(ceiling(k)==floor(k))
        {
            #### Determine the acceptance probabilities
            accept.lambda <- 100 * accept[1] / accept[2]
            accept.delta <- 100 * accept[3] / accept[4]
            accept.phi <- 100 * accept[5] / accept[6]            
            if(!is.null(X))
            {
                accept.beta <- 100 * accept[7] / accept[8]
                if(accept.beta > 40)
                {
                    proposal.sd.beta <- proposal.sd.beta + 0.1 * proposal.sd.beta
                }else if(accept.beta < 20)              
                {
                    proposal.sd.beta <- proposal.sd.beta - 0.1 * proposal.sd.beta
                }else
                {
                }   
                accept.all <- accept.all + accept
                accept <- rep(0,8)
            }else
            {
                accept.all <- accept.all + accept
                accept <- rep(0,6)   
            } 
            
            
            
            
            #### lambda tuning parameter
            if(accept.lambda > 40)
            {
                proposal.sd.lambda <- min(proposal.sd.lambda + 0.1 * proposal.sd.lambda, 10)
            }else if(accept.lambda < 20)              
            {
                proposal.sd.lambda <- proposal.sd.lambda - 0.1 * proposal.sd.lambda
            }else
            {
            }
            
            
            #### delta tuning parameter               
            if(accept.delta > 50)
            {
                proposal.sd.delta <- min(proposal.sd.delta + 0.1 * proposal.sd.delta, 10)
            }else if(accept.delta < 40)              
            {
                proposal.sd.delta <- proposal.sd.delta - 0.1 * proposal.sd.delta
            }else
            {
            }         
            
            
            #### phi tuning parameter
            if(accept.phi > 50)
            {
                proposal.sd.phi <- proposal.sd.phi + 0.1 * proposal.sd.phi
            }else if(accept.phi < 40)              
            {
                proposal.sd.phi <- proposal.sd.phi - 0.1 * proposal.sd.phi
            }else
            {
            }
        }else
        {   
        }
        
        
        
        ################################       
        ## print progress to the console
        ################################
        if(j %in% percentage.points & verbose)
        {
            setTxtProgressBar(progressBar, j/n.sample)
        }
    }
    
    # end timer
    if(verbose)
    {
        cat("\nSummarising results")
        close(progressBar)
    }else
    {}
    
    
    ###################################
    #### Summarise and save the results 
    ###################################
    ## Compute the acceptance rates
    accept.lambda <- 100 * accept.all[1] / accept.all[2]
    accept.delta <- 100 * accept.all[3] / accept.all[4]
    accept.phi <- 100 * accept.all[5] / accept.all[6]
    accept.gamma <- 100
    if(!is.null(X))
    {
        accept.beta <- 100 * accept.all[7] / accept.all[8]   
        accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.gamma)
        names(accept.final) <- c("beta", "lambda", "delta", "phi", "rho.T")   
    }else
    {
        accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.gamma)
        names(accept.final) <- c("lambda", "delta", "phi", "rho.T")   
    }
    
    
    
    ## DIC
    median.Z <- round(apply(samples.Z,2,median), 0)       
    median.lambda <- apply(samples.lambda, 2, median)
    median.mu <- matrix(median.lambda[median.Z], nrow=K, ncol=N, byrow=FALSE)
    
    if(!is.null(X))
    {
        median.beta <- apply(samples.beta,2,median)
        regression.mat <- matrix(X.standardised %*% median.beta, nrow=K, ncol=N, byrow=FALSE)     
    }else
    {
    }
    
    median.phi <- matrix(apply(samples.phi, 2, median), nrow=K, byrow=FALSE)
    lp.median <- as.numeric(median.mu + offset.mat + median.phi + regression.mat)   
    median.prob <- exp(lp.median)  / (1 + exp(lp.median))
    fitted.median <- trials * median.prob
    deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE))
    p.d <- median(samples.deviance) - deviance.fitted
    DIC <- 2 * median(samples.deviance) - deviance.fitted     
    
    
    #### Watanabe-Akaike Information Criterion (WAIC)
    LPPD <- sum(log(apply(samples.like,2,mean)), na.rm=TRUE)
    p.w <- sum(apply(log(samples.like),2,var), na.rm=TRUE)
    WAIC <- -2 * (LPPD - p.w)
    
    
    ## Compute the LMPL
    CPO <- rep(NA, N.all)
    for(j in 1:N.all)
    {
        CPO[j] <- 1/median((1 / dbinom(x=Y[j], size=trials[j], prob=(samples.fitted[ ,j] / trials[j]))))    
    }
    LMPL <- sum(log(CPO))  
    
    
    ## Create the Fitted values
    fitted.values <- apply(samples.fitted, 2, median)
    residuals <- as.numeric(Y) - fitted.values
    
    
    #### transform the parameters back to the origianl covariate scale.
    if(!is.null(X))
    {    
        samples.beta.orig <- samples.beta
        number.cts <- sum(X.indicator==1)     
        if(number.cts>0)
        {
            for(r in 1:p)
            {
                if(X.indicator[r]==1)
                {
                    samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
                }else
                {
                }
            }
        }else
        {
        }
    }else
    {}
    
    
    
    #### Create a summary object
    summary.hyper <- array(NA, c(3, 7))     
    summary.hyper[1,1:3] <- quantile(samples.delta, c(0.5, 0.025, 0.975))
    summary.hyper[2,1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
    summary.hyper[3,1:3] <- quantile(samples.gamma, c(0.5, 0.025, 0.975))
    rownames(summary.hyper) <- c("delta", "tau2", "rho.T")      
    summary.hyper[1, 4:7] <- c(n.keep, accept.delta, effectiveSize(mcmc(samples.delta)), geweke.diag(mcmc(samples.delta))$z)   
    summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2)), geweke.diag(mcmc(samples.tau2))$z)   
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.gamma)), geweke.diag(mcmc(samples.gamma))$z)   
    
    summary.lambda <- array(NA, c(G,1))
    summary.lambda <- t(apply(samples.lambda, 2, quantile, c(0.5, 0.025, 0.975)))
    summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(mcmc(samples.lambda)), geweke.diag(mcmc(samples.lambda))$z)
    summary.lambda <- matrix(summary.lambda, ncol=7)
    rownames(summary.lambda) <- paste("lambda", 1:G, sep="")
    
    
    if(!is.null(X))
    {
        samples.beta.orig <- mcmc(samples.beta.orig)
        summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
        summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
        rownames(summary.beta) <- colnames(X)
        colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
        summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)    
    }else
    {
        summary.results <- rbind(summary.lambda, summary.hyper)    
    }
    
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    colnames(summary.results) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")    
    
    
    
    ## Compile and return the results
    modelfit <- c(DIC, p.d, WAIC, p.w, LMPL)
    names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL")
    if(is.null(X)) samples.beta.orig = NA

    samples <- list(beta=mcmc(samples.beta.orig), lambda=mcmc(samples.lambda),  Z=mcmc(samples.Z), delta=mcmc(samples.delta), phi = mcmc(samples.phi), tau2=mcmc(samples.tau2), rho.T=mcmc(samples.gamma), fitted=mcmc(samples.fitted), deviance=mcmc(samples.deviance))
    model.string <- c("Likelihood model - Binomial (logit link function)", "\nLatent structure model - Localised autoregressive CAR model\n")
    results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=median.Z, formula=formula, model=model.string,  X=X)
    class(results) <- "carbayesST"
    if(verbose)
    {
        b<-proc.time()
        cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
    return(results)
}
