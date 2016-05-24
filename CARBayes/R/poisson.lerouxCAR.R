poisson.lerouxCAR <- function(formula, data=NULL,  W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, fix.rho=FALSE, rho=NULL, verbose=TRUE)
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
    
    n <- nrow(X)
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
    
    
    
    #### Response variable
    ## Create the response
    Y <- model.response(frame)
    which.miss <- as.numeric(!is.na(Y))
    n.miss <- n - sum(which.miss)
    Y.miss <- Y
    Y.miss[which.miss==0] <- median(Y, na.rm=TRUE)
    
    ## Check for errors
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
    int.check <- n - n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    
    
    #### Offset variable
    ## Create the offset
    offset <- try(model.offset(frame), silent=TRUE)
    
    ## Check for errors
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,n)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    ## Check for errors on rho and fix.rho
    if(!is.logical(fix.rho)) stop("fix.rho is not logical.", call.=FALSE)   
    if(fix.rho & is.null(rho)) stop("rho is fixed but an initial value was not set.", call.=FALSE)   
    if(fix.rho & !is.numeric(rho) ) stop("rho is fixed but is not numeric.", call.=FALSE)  
    
    #### Initial parameter values
    mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    
    log.Y <- log(Y)
    log.Y[Y==0] <- -0.1  
    res.temp <- log.Y - X.standardised %*% beta.mean - offset
    res.sd <- sd(res.temp, na.rm=TRUE)/5
    phi <- rnorm(n=n, mean=rep(0,n), sd=res.sd)
    tau2 <- var(phi) / 10
    
    if(!fix.rho) rho <- runif(1)
    if(rho<0 ) stop("rho is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho is outside the range [0, 1].", call.=FALSE) 
    
    
    #### Priors
    ## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
    
    
    ## Checks    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
    
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)
    
    if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
    
    
    
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
    samples.phi <- array(NA, c(n.keep, n))
    samples.tau2 <- array(NA, c(n.keep, 1))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
    samples.deviance <- array(NA, c(n.keep, 1))
    samples.like <- array(NA, c(n.keep, n))
    samples.fitted <- array(NA, c(n.keep, n))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    ## Metropolis quantities
    accept.all <- rep(0,6)
    accept <- accept.all
    proposal.sd.beta <- 0.01
    proposal.sd.phi <- 0.1
    proposal.sd.rho <- 0.02
    proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
    chol.proposal.corr.beta <- chol(proposal.corr.beta) 
    tau2.posterior.shape <- prior.tau2[1] + 0.5 * n
    
    
    
    #### CAR quantities
    ## Checks
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    
    if(fix.rho & rho==0)
    {
        ## Set up a dummy W matrix to use in the code as it will not affect the results
        W <- array(0, c(n,n))
        for(r in 2:n)
        {
            W[(r-1), r] <- 1   
            W[r, (r-1)] <- 1
        }
    }else
    {
        if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
    }
    
    
    ## Create the triplet object
    W.triplet <- c(NA, NA, NA)
    for(i in 1:n)
    {
        for(j in 1:n)
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
    
    
    ## Create the start and finish points for W updating
    W.begfin <- array(NA, c(n, 2))     
    temp <- 1
    for(i in 1:n)
    {
        W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
        temp <- temp + n.neighbours[i]
    }
    
    
    ## Create the determinant     
    if(!fix.rho)
    {
        Wstar <- diag(apply(W,1,sum)) - W
        Wstar.eigen <- eigen(Wstar)
        Wstar.val <- Wstar.eigen$values
        det.Q <- 0.5 * sum(log((rho * Wstar.val + (1-rho))))    
    }else
    {}   
    
    
    ## Check for islands
    W.list<- mat2listw(W)
    W.nb <- W.list$neighbours
    W.islands <- n.comp.nb(W.nb)
    islands <- W.islands$comp.id
    n.islands <- max(W.islands$nc)
    if(rho==1) tau2.posterior.shape <- prior.tau2[1] + 0.5 * (n-n.islands)   

        
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
            prob <- poissonbetaupdate(X.standardised, n, p, beta, proposal.beta, offset.temp, Y.miss, prior.mean.beta, prior.var.beta, which.miss)
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
        
        
        
        ####################
        ## Sample from phi
        ####################
        beta.offset <- X.standardised %*% beta + offset
        temp1 <- poissoncarupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=n, phi=phi, tau2=tau2, y=Y.miss, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset, which.miss)
        phi <- temp1[[1]]
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }
        else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
        accept[3] <- accept[3] + temp1[[2]]
        accept[4] <- accept[4] + n    
        
        
        
        ##################
        ## Sample from tau2
        ##################
        temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, n, phi, phi, rho)
        tau2.posterior.scale <- temp2 + prior.tau2[2] 
        tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
        
        
        
        ##################
        ## Sample from rho
        ##################
        if(!fix.rho)
        {
        proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)   
        temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, n, phi, phi, proposal.rho)
        det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
        logprob.current <- det.Q - temp2 / tau2
        logprob.proposal <- det.Q.proposal - temp3 / tau2
        prob <- exp(logprob.proposal - logprob.current)
        
        #### Accept or reject the proposal
        if(prob > runif(1))
        {
            rho <- proposal.rho
            det.Q <- det.Q.proposal
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
        fitted <- exp(lp)
        deviance.all <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)
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
            samples.tau2[ele, ] <- tau2
            if(!fix.rho) samples.rho[ele, ] <- rho
            samples.deviance[ele, ] <- deviance
            samples.like[ele, ] <- like
            samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- rpois(n=n.miss, lambda=fitted[which.miss==0])
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
    accept.tau2 <- 100
    accept.final <- c(accept.beta, accept.phi, accept.rho, accept.tau2)
    names(accept.final) <- c("beta", "phi", "rho", "tau2")
    
    
    ## Deviance information criterion (DIC)
    median.beta <- apply(samples.beta, 2, median)
    median.phi <- apply(samples.phi, 2, median)
    fitted.median <- exp(X.standardised %*% median.beta + median.phi + offset)
    deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.median, log=TRUE), na.rm=TRUE)
    p.d <- median(samples.deviance) - deviance.fitted
    DIC <- 2 * median(samples.deviance) - deviance.fitted     
    
    
    #### Watanabe-Akaike Information Criterion (WAIC)
    LPPD <- sum(log(apply(samples.like,2,mean)), na.rm=TRUE)
    p.w <- sum(apply(log(samples.like),2,var), na.rm=TRUE)
    WAIC <- -2 * (LPPD - p.w)
    
    #### Compute the Conditional Predictive Ordinate
    CPO <- rep(NA, n)
    for(j in 1:n)
    {
        CPO[j] <- 1/median((1 / dpois(x=Y[j], lambda=samples.fitted[ ,j])))    
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
    
    
    summary.hyper <- array(NA, c(2 ,7))
    summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
    summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)
    if(!fix.rho)
    {
        summary.hyper[2, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
        summary.hyper[2, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
        summary.hyper[2, 1:3] <- c(rho, rho, rho)
        summary.hyper[2, 4:7] <- rep(NA, 4)
    }
    
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "rho")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    
    
    #### Create the Fitted values and residuals
    fitted.values <- apply(samples.fitted, 2, median)
    residuals <- as.numeric(Y) - fitted.values
    
    
    ## Compile and return the results
    modelfit <- c(DIC, p.d, WAIC, p.w, LMPL)
    names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL")
    model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Leroux CAR\n")
    if(fix.rho) samples.rho=NA
    if(n.miss==0) samples.Y = NA
    
    samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), rho=mcmc(samples.rho), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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
