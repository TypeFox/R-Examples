binomial.MVlerouxCAR <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, fix.rho=FALSE, rho=NULL, verbose=TRUE)
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
    if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
    
    
    
    #### Design matrix
    ## Create the matrix
    X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
    if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
    if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
    
    N.all <- nrow(X)
    p <- ncol(X)
    
    ## Check for linearly related columns
    cor.X <- suppressWarnings(cor(X))
    diag(cor.X) <- 0
    
    if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    
    if(p>1)
    {
        if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
    }else
    {
    }
    
    ## Standardise the matrix
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
        }else if(length(table(X[ ,j]))==1)
        {
            X.indicator[j] <- 2
        }else
        {
            X.indicator[j] <- 0
        }
    }
    
    
    
    #### Response variable and trials
    ## Create the response
    Y <- model.response(frame)
    failures <- trials - Y
    
    ## Check for errors
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
    int.check <- N.all-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
    
    
    ## Create the missing value indicator
    which.miss <- as.numeric(!is.na(Y))
    n.miss <- N.all - sum(which.miss)
    Y.miss <- Y
    Y.miss[which.miss==0] <- median(Y, na.rm=TRUE)
    failures.miss <- failures
    failures.miss[which.miss==0] <- median(failures, na.rm=TRUE)
    
    
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
    int.check <- N.all - n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)
    
    
    
    #### Offset variable
    ## Create the offset
    offset <- try(model.offset(frame), silent=TRUE)
    
    ## Check for errors
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,N.all)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    #### W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    K <- nrow(W)
    if(ceiling(N.all/K)!= floor(N.all/K)) stop("The number of data points divided by the number of rows in W is not a whole number.", call.=FALSE)
    J <- N.all / K
    if(ncol(W)!= K) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(!is.symmetric.matrix(W)) stop("W is not symmetric.", call.=FALSE)


    
    ## Check for errors on rho and fix.rho
    if(!is.logical(fix.rho)) stop("fix.rho is not logical.", call.=FALSE)   
    if(fix.rho & is.null(rho)) stop("rho is fixed but an initial value was not set.", call.=FALSE)   
    if(fix.rho & !is.numeric(rho) ) stop("rho is fixed but is not numeric.", call.=FALSE)  
    
    #### Initial parameter values
    ## Regression parameters beta
    dat <- cbind(Y, failures)
    mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
    
    theta.hat <- Y / trials
    theta.hat[theta.hat==0] <- 0.01
    theta.hat[theta.hat==1] <- 0.99
    res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta.mean - offset
    res.sd <- sd(res.temp, na.rm=TRUE)/5
    phi <- rnorm(n=N.all, mean=rep(0,N.all), sd=res.sd)
    phi.mat <- matrix(phi, nrow=K, byrow=TRUE)
    Sigma <- cov(phi.mat)
    Sigma.inv <- solve(Sigma)
    if(!fix.rho) rho <- runif(1)
    if(rho<0 ) stop("rho is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho is outside the range [0, 1].", call.=FALSE)  
    
    
    if(fix.rho & rho==0)
    {
        ## Set up a dummy W matrix to use in the code as it will not affect the results
        W <- array(0, c(K,K))
        for(r in 2:K)
        {
            W[(r-1), r] <- 1   
            W[r, (r-1)] <- 1
        }
    }else
    {
        if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
    }
    
    
    #### Priors
    ## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.Sigma.df)) prior.Sigma.df <- J+1
    if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- diag(rep(1,J))
    
    ## Checks    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
    
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)
    
    if(length(prior.Sigma.df)!=1) stop("the prior value for prior.Sigma.df is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.Sigma.df)) stop("the prior value for prior.Sigma.df is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.Sigma.df))!=0) stop("the prior value for prior.Sigma.df has missing values.", call.=FALSE)    
    
    if(nrow(prior.Sigma.scale)!=J) stop("prior.Sigma.scale is the wrong dimension.", call.=FALSE)    
    if(ncol(prior.Sigma.scale)!=J) stop("prior.Sigma.scale is the wrong dimension.", call.=FALSE)    
    if(!is.numeric(prior.Sigma.scale)) stop("prior.Sigma.scale has non-numeric values.", call.=FALSE)    
    if(sum(is.na(prior.Sigma.scale))!=0) stop("prior.Sigma.scale has missing values.", call.=FALSE)    
    if(!is.positive.definite(prior.Sigma.scale)) stop("prior.Sigma.scale is not a positive definite matrix.", call.=FALSE)
    if(!is.symmetric.matrix(prior.Sigma.scale)) stop("prior.Sigma.scale is not symmetric.", call.=FALSE)
    
    
    #### MCMC quantities
    ## Checks
    if(is.null(burnin)) stop("the burnin argument is missing", call.=FALSE)
    if(is.null(n.sample)) stop("the n.sample argument is missing", call.=FALSE)
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(n.sample <= thin)  stop("thin is greater than n.sample.", call.=FALSE)
    if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
    if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
    if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 
    
    
    ## Compute the blocking structure for beta     
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
    
    
    
   

    
    ## Matrices to store samples
    n.keep <- floor((n.sample - burnin)/thin)
    samples.beta <- array(NA, c(n.keep, p))
    samples.phi <- array(NA, c(n.keep, N.all))
    samples.Sigma <- array(NA, c(n.keep, J, J))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
    samples.deviance <- array(NA, c(n.keep, 1))
    samples.like <- array(NA, c(n.keep, N.all))
    samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
    ## Metropolis quantities
    accept.all <- rep(0,6)
    accept <- accept.all
    proposal.sd.beta <- 0.01
    proposal.sd.phi <- 0.1
    proposal.sd.rho <- 0.02
    proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
    chol.proposal.corr.beta <- chol(proposal.corr.beta) 
    Sigma.post.df <- prior.Sigma.df + K   
    
    #### CAR quantities
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
    n.triplet <- nrow(W.triplet) 
    W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
    n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
    Wstar <- diag(apply(W,1,sum)) - W
    Q <- rho * Wstar + diag(rep(1-rho,K))
    
    ## Create the start and finish points for W updating
    W.begfin <- array(NA, c(K, 2))     
    temp <- 1
    for(i in 1:K)
    {
        W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
        temp <- temp + n.neighbours[i]
    }
    
    
    ## Create the determinant     
        if(!fix.rho)
        {
        Wstar.eigen <- eigen(Wstar)
        Wstar.val <- Wstar.eigen$values
        det.Q <- sum(log((rho * Wstar.val + (1-rho))))    
        }else
        {}
    
    ## Check for islands
    W.list<- mat2listw(W)
    W.nb <- W.list$neighbours
    W.islands <- n.comp.nb(W.nb)
    islands <- W.islands$comp.id
    islands.all <- rep(islands,J)
    n.islands <- max(W.islands$nc)
    if(rho==1) Sigma.post.df <- prior.Sigma.df + K - n.islands   

     
    ## Make matrix versions of the variables
    offset.mat <- matrix(offset, nrow=K, ncol=J, byrow=TRUE) 
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=J, byrow=TRUE)
    Y.mat <- matrix(Y, nrow=K, ncol=J, byrow=TRUE)
    failures.mat <- matrix(failures, nrow=K, ncol=J, byrow=TRUE)
    Y.mat.miss <- matrix(Y.miss, nrow=K, ncol=J, byrow=TRUE)
    failures.mat.miss <- matrix(failures.miss, nrow=K, ncol=J, byrow=TRUE)
    which.miss.mat <- matrix(which.miss, nrow=K, ncol=J, byrow=TRUE)
    
    
    ###########################
    #### Run the Bayesian model
    ###########################
    ## Start timer
    if(verbose)
    {
        cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n", sep = " ")
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
        proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
        proposal.beta <- beta
        offset.temp <- phi + offset
        
        for(r in 1:n.beta.block)
        {
            proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
            prob <- binomialbetaupdate(X.standardised, N.all, p, beta, proposal.beta, offset.temp, Y.miss, failures.miss, prior.mean.beta, prior.var.beta, which.miss)
            if(prob > runif(1))
            {
                beta[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
                accept[1] <- accept[1] + 1  
            }else
            {
                proposal.beta[beta.beg[r]:beta.fin[r]] <- beta[beta.beg[r]:beta.fin[r]]
            }
        }
        
        accept[2] <- accept[2] + n.beta.block    
        regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=J, byrow=TRUE)   
        
        
        ####################
        ## Sample from phi
        ####################
        chol.Sigma <- t(chol(proposal.sd.phi * Sigma))
        rand <- matrix(rnorm(n=N.all), nrow=K)
        den.offset <- rho * W.triplet.sum + 1 - rho
        phi.offset <- regression.mat + offset.mat
        temp1 <- binomialmcarupdate(W.triplet, W.begfin, W.triplet.sum,  K, J, phi.mat, Y.mat.miss, failures.mat.miss, phi.offset, den.offset, Sigma, Sigma.inv, rho,  chol.Sigma, rand, which.miss.mat)      
            if(rho<1)
            {
            phi.mat <- temp1[[1]] - mean(temp1[[1]])
            phi <- as.numeric(t(phi.mat))
            }else
            {
            phi.mat <- temp1[[1]]
            phi <- as.numeric(t(phi.mat))
            phi[which(islands.all==1)] <- phi[which(islands.all==1)] - mean(phi[which(islands.all==1)]) 
            phi.mat <- matrix(phi, nrow=K, byrow=TRUE)
            }
        accept[3] <- accept[3] + temp1[[2]]
        accept[4] <- accept[4] + K    
        
        
   
        
        ##################
        ## Sample from Sigma
        ##################
        Sigma.post.scale <- t(phi.mat) %*% Q %*% phi.mat + prior.Sigma.scale
        Sigma <- riwish(Sigma.post.df, Sigma.post.scale)
        Sigma.inv <- solve(Sigma)
        
        
        ##################
        ## Sample from rho
        ##################
        if(!fix.rho)
        {
        ## Propose a new value
        proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)
        Q.prop <- proposal.rho * Wstar + diag(rep(1-proposal.rho), K)
        det.Q.prop <-  sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))    
            
        ## Compute the acceptance rate
        logprob.current <- 0.5 * J * det.Q - 0.5 * sum(diag(t(phi.mat) %*% Q %*% phi.mat %*% Sigma.inv))
        logprob.proposal <- 0.5 * J * det.Q.prop - 0.5 * sum(diag(t(phi.mat) %*% Q.prop %*% phi.mat %*% Sigma.inv))
        prob <- exp(logprob.proposal - logprob.current)
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q <- det.Q.prop
            Q <- Q.prop
            accept[5] <- accept[5] + 1           
            }else
            {
            }              
            accept[6] <- accept[6] + 1       
        }else
        {}
        
 
        
        #########################
        ## Calculate the deviance
        #########################
        lp <- as.numeric(X.standardised %*% beta) + phi + offset
        prob <- exp(lp)  / (1 + exp(lp))
        fitted <- trials * prob
        deviance.all <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)
        like <- exp(deviance.all)
        deviance <- -2 * sum(deviance.all, na.rm=TRUE)
        
        
        
        ###################
        ## Save the results
        ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
            ele <- (j - burnin) / thin
            samples.beta[ele, ] <- beta
            samples.phi[ele, ] <- phi
            samples.Sigma[ele, , ] <- Sigma
            if(!fix.rho) samples.rho[ele, ] <- rho
            samples.deviance[ele, ] <- deviance
            samples.like[ele, ] <- like
            samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- rbinom(n=n.miss, size=trials[which.miss==0], prob=prob[which.miss==0])
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
            accept.beta <- 100 * accept[1] / accept[2]
            accept.phi <- 100 * accept[3] / accept[4]
            accept.rho <- 100 * accept[5] / accept[6]
            if(is.na(accept.rho)) accept.rho <- 45
            accept.all <- accept.all + accept
            accept <- c(0,0,0,0,0,0)
            
            #### beta tuning parameter
            if(accept.beta > 40)
            {
                proposal.sd.beta <- proposal.sd.beta + 0.1 * proposal.sd.beta
            }else if(accept.beta < 20)              
            {
                proposal.sd.beta <- proposal.sd.beta - 0.1 * proposal.sd.beta
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
            
            #### rho tuning parameter
            if(accept.rho > 50)
            {
                proposal.sd.rho <- min(proposal.sd.rho + 0.1 * proposal.sd.rho, 0.5)
            }else if(accept.rho < 40)              
            {
                proposal.sd.rho <- proposal.sd.rho - 0.1 * proposal.sd.rho
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
    accept.beta <- 100 * accept.all[1] / accept.all[2]
    accept.phi <- 100 * accept.all[3] / accept.all[4]
    if(!fix.rho)
    {
        accept.rho <- 100 * accept.all[5] / accept.all[6]
    }else
    {
        accept.rho <- NA    
    }
    accept.Sigma <- 100
    accept.final <- c(accept.beta, accept.phi, accept.rho, accept.Sigma)
    names(accept.final) <- c("beta", "phi", "rho", "Sigma")
    
    ## Deviance information criterion (DIC)
    median.beta <- apply(samples.beta, 2, median)
    median.phi <- apply(samples.phi, 2, median)
    median.logit <- as.numeric(X.standardised %*% median.beta) + median.phi + offset    
    median.prob <- exp(median.logit)  / (1 + exp(median.logit))
    fitted.median <- trials * median.prob
    deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE), na.rm=TRUE)
    p.d <- median(samples.deviance) - deviance.fitted
    DIC <- 2 * median(samples.deviance) - deviance.fitted     
    
    
    #### Watanabe-Akaike Information Criterion (WAIC)
    LPPD <- sum(log(apply(samples.like,2,mean)), na.rm=TRUE)
    p.w <- sum(apply(log(samples.like),2,var), na.rm=TRUE)
    WAIC <- -2 * (LPPD - p.w)
    
    
    #### Compute the Conditional Predictive Ordinate
    CPO <- rep(NA, N.all)
    for(j in 1:N.all)
    {
        CPO[j] <- 1/median((1 / dbinom(x=Y[j], size=trials[j], prob=(samples.fitted[ ,j] / trials[j]))))    
    }
    LMPL <- sum(log(CPO), na.rm=TRUE)     
    
    
    
    
    #### transform the parameters back to the origianl covariate scale.
    samples.beta.orig <- samples.beta
    number.cts <- sum(X.indicator==1)     
    if(number.cts>0)
    {
        for(r in 1:p)
        {
            if(X.indicator[r]==1)
            {
                samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
            }else if(X.indicator[r]==2 & p>1)
            {
                X.transformed <- which(X.indicator==1)
                samples.temp <- as.matrix(samples.beta[ ,X.transformed])
                for(s in 1:length(X.transformed))
                {
                    samples.temp[ ,s] <- samples.temp[ ,s] * X.mean[X.transformed[s]]  / X.sd[X.transformed[s]]
                }
                intercept.adjustment <- apply(samples.temp, 1,sum) 
                samples.beta.orig[ ,r] <- samples.beta[ ,r] - intercept.adjustment
            }else
            {
            }
        }
    }else
    {
    }
    
    
    
    #### Create a summary object
    samples.beta.orig <- mcmc(samples.beta.orig)
    summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
    summary.hyper <- array(NA, c((J+1) ,7))
    summary.hyper[1:J, 1] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.5)))
    summary.hyper[1:J, 2] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.025)))
    summary.hyper[1:J, 3] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.975)))
    summary.hyper[1:J, 4] <- n.keep
    summary.hyper[1:J, 5] <- accept.Sigma
    summary.hyper[1:J, 6] <- diag(apply(samples.Sigma, c(2,3), effectiveSize))
        for(r in 1:J)
        {
        summary.hyper[r, 7] <- geweke.diag(samples.Sigma[ ,r,r])$z    
        }
    
    if(!fix.rho)
    {
        summary.hyper[(J+1), 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
        summary.hyper[(J+1), 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
        summary.hyper[(J+1), 1:3] <- c(rho, rho, rho)
        summary.hyper[(J+1), 4:7] <- rep(NA, 4)
    }
    
    
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[(p+1): nrow(summary.results)] <- c(paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    
    
    
    #### Create the Fitted values and residuals
    fitted.values <- apply(samples.fitted, 2, median)
    residuals <- as.numeric(Y) - fitted.values
    
    
    ## Compile and return the results
    modelfit <- c(DIC, p.d, WAIC, p.w, LMPL)
    names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL")
    model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Leroux MCAR\n")
    if(fix.rho) samples.rho=NA
    if(n.miss==0) samples.Y = NA
    samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), Sigma=samples.Sigma, rho=mcmc(samples.rho), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
    results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, X=X)
    class(results) <- "carbayes"
    
    if(verbose)
    {
        b<-proc.time()
        cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
    return(results)
}

