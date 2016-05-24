binomial.CARar <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, fix.rho.S=FALSE, rho.S=NULL, fix.rho.T=FALSE, rho.T=NULL, verbose=TRUE)
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
    
    
    #### Format and check the neighbourhood matrix W
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    K <- nrow(W)
    if(ncol(W)!= K) stop("W has the wrong number of columns.", call.=FALSE)
    if(nrow(W)!= K) stop("W has the wrong number of rows.", call.=FALSE)
    if(floor(N.all/K)!=ceiling(N.all/K)) stop("The number of spatial areas is not a multiple of the number of data points.", call.=FALSE)
    N <- N.all / K
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)   
    
    #### Response variable
Y <- model.response(frame)
Y <- as.numeric(Y)
failures <- trials - Y

## Check for errors
if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- N.all-sum(ceiling(trials)==floor(trials))
if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)

## Create the missing value indicator
which.miss <- as.numeric(!is.na(Y))
which.miss.mat <- matrix(which.miss, nrow=K, ncol=N, byrow=FALSE)
n.miss <- N.all - sum(which.miss)
Y.miss <- Y
Y.miss[which.miss==0] <- median(Y, na.rm=TRUE)
failures.miss <- failures
failures.miss[which.miss==0] <- median(failures, na.rm=TRUE)

if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
int.check <- N.all- n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
if(int.check > 0) stop("the response variable has non-integer values.", call.=FALSE)
if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


    #### Offset variable
    offset <- try(model.offset(frame), silent=TRUE)
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,N.all)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    
 

#### Specify the initial parameter values
dat <- cbind(Y, failures)
mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    
theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=N.all, mean=0, sd = res.sd)
tau2 <- var(phi)/10
    if(fix.rho.S)
    {
        rho <- rho.S
    }else
    {
        rho <- runif(1)       
    }
    
    if(fix.rho.T)
    {
        gamma <- rho.T
    }else
    {
        gamma <- runif(1)       
    }   
    
    
    #### Check and specify the priors
    ## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
    
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
    
    
    ## Check for errors on rho and fix.rho
    if(!is.logical(fix.rho.S)) stop("fix.rho.S is not logical.", call.=FALSE)   
    if(fix.rho.S & is.null(rho.S)) stop("rho.S is fixed but an initial value was not set.", call.=FALSE)   
    if(fix.rho.S & !is.numeric(rho.S) ) stop("rho.S is not numeric.", call.=FALSE)  
    if(rho<0 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)  
    
    ## Check for errors on rho and fix.rho
    if(!is.logical(fix.rho.T)) stop("fix.rho.T is not logical.", call.=FALSE)   
    if(fix.rho.T & is.null(rho.T)) stop("rho.T is fixed but an initial value was not set.", call.=FALSE)   
    if(fix.rho.T & !is.numeric(rho.T) ) stop("rho.T is not numeric.", call.=FALSE)  
    if(gamma<0 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
    if(gamma>1 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
    
    #### Set up matrices to store samples
    n.keep <- floor((n.sample - burnin)/thin)
    samples.beta <- array(NA, c(n.keep, p))
    samples.phi <- array(NA, c(n.keep, N.all))
    samples.tau2 <- array(NA, c(n.keep, 1))
    if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.T) samples.gamma <- array(NA, c(n.keep, 1))
    samples.fitted <- array(NA, c(n.keep, N.all))
    samples.like <- array(NA, c(n.keep, N.all))
    samples.deviance <- array(NA, c(n.keep, 1))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
    #### Specify the Metropolis quantities
    accept.all <- rep(0,6)
    accept <- accept.all
    proposal.sd.phi <- 0.1
    proposal.sd.rho <- 0.05
    proposal.sd.beta <- 0.01
    proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
    chol.proposal.corr.beta <- chol(proposal.corr.beta)     
    tau2.shape <- prior.tau2[1] + N.all/2
    

    
    
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
    
    
    ## Create the determinant     
    if(!fix.rho.S) 
    {
        Wstar <- diag(apply(W,1,sum)) - W
        Wstar.eigen <- eigen(Wstar)
        Wstar.val <- Wstar.eigen$values
        det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
    }else
    {}


     
#### Specify quantities that do not change
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE)
Y.mat.miss <- matrix(Y.miss, nrow=K, ncol=N, byrow=FALSE)
failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=FALSE)
failures.mat.miss <- matrix(failures.miss, nrow=K, ncol=N, byrow=FALSE)
phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)   

 
## Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
if(rho==1 & gamma==1) 
{
    tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * (K-1))/2
}else if(rho==1)
{
    tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + (N * (K-1))/2        
}else if(gamma==1)
{
    tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * K)/2          
}else
{}



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
        proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
        proposal.beta <- beta
        offset.temp <- as.numeric(offset.mat + phi.mat)     
        
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
        regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
        

        
        ####################
        ## Sample from phi
        ####################
        phi.offset <- offset.mat + regression.mat
        den.offset <- rho * W.triplet.sum + 1 - rho
        temp1 <- binomialarcarupdate(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, gamma, rho, Y.mat.miss, failures.mat.miss, proposal.sd.phi, phi.offset, den.offset, which.miss.mat)      
        phi.temp <- temp1[[1]]
        phi <- as.numeric(phi.temp)  - mean(as.numeric(phi.temp))
        phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
        accept[3] <- accept[3] + temp1[[2]]
        accept[4] <- accept[4] + K*N    
        
        

        ####################
        ## Sample from gamma
        ####################
        if(!fix.rho.T)
        {
        temp2 <- gammaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho)
        mean.gamma <- temp2[[1]] / temp2[[2]]
        sd.gamma <- sqrt(tau2 / temp2[[2]])
        gamma <- rtrunc(n=1, spec="norm", a=0, b=1, mean=mean.gamma, sd=sd.gamma)  
        }else
        {}
        

                
        ####################
        ## Samples from tau2
        ####################
        temp3 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho, gamma)
        tau2.scale <- temp3 + prior.tau2[2] 
        tau2 <- 1 / rgamma(1, tau2.shape, scale=(1/tau2.scale)) 
        
        

        ##################
        ## Sample from rho
        ##################
        if(!fix.rho.S)
        {
        proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)   
        temp4 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, proposal.rho, gamma)
        det.Q.W.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
        logprob.current <- N * det.Q.W - temp3 / tau2
        logprob.proposal <- N * det.Q.W.proposal - temp4 / tau2
        prob <- exp(logprob.proposal - logprob.current)
        if(prob > runif(1))
        {
            rho <- proposal.rho
            det.Q.W <- det.Q.W.proposal
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
        lp <- as.numeric(offset.mat + regression.mat + phi.mat)
        prob <- exp(lp) / (1+exp(lp))
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
            samples.phi[ele, ] <- as.numeric(phi)
            if(!fix.rho.S) samples.rho[ele, ] <- rho
            if(!fix.rho.T) samples.gamma[ele, ] <- gamma
            samples.tau2[ele, ] <- tau2
            samples.deviance[ele, ] <- deviance
            samples.fitted[ele, ] <- fitted
            samples.like[ele, ] <- like
            if(n.miss>0) samples.Y[ele, ] <- rbinom(n=n.miss, size=trials[which.miss==0], prob=prob[which.miss==0])
        }else{
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
            accept <- rep(0,6)
            
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
if(!fix.rho.S)
{
    accept.rho <- 100 * accept.all[5] / accept.all[6]
}else
{
    accept.rho <- NA    
}
accept.final <- c(accept.beta, accept.phi, accept.rho, 100)
names(accept.final) <- c("beta", "phi", "rho.S", "rho.T")

    
    
## Compute information criterion (DIC, DIC3, WAIC)
median.beta <- apply(samples.beta,2,median)
regression.mat <- matrix(X.standardised %*% median.beta, nrow=K, ncol=N, byrow=FALSE)   
median.phi <- matrix(apply(samples.phi, 2, median), nrow=K, ncol=N)
lp.median <- as.numeric(offset.mat + median.phi + regression.mat)   
median.prob <- exp(lp.median)  / (1 + exp(lp.median))
fitted.median <- trials * median.prob
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE), na.rm=TRUE)
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
LMPL <- sum(log(CPO), na.rm=TRUE)  
    

## Create the Fitted values
fitted.values <- apply(samples.fitted, 2, median)
residuals <- as.numeric(Y) - fitted.values
    
    
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
    
summary.hyper <- array(NA, c(3, 7))    
rownames(summary.hyper) <- c("tau2", "rho.S", "rho.T")     
summary.hyper[1,1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2)), geweke.diag(mcmc(samples.tau2))$z)     

if(!fix.rho.S)
{
    summary.hyper[2, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[2, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
}else
{
    summary.hyper[2, 1:3] <- c(rho, rho, rho)
    summary.hyper[2, 4:7] <- rep(NA, 4)
}

if(!fix.rho.T)
{
    summary.hyper[3, 1:3] <- quantile(samples.gamma, c(0.5, 0.025, 0.975))
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.gamma)), geweke.diag(mcmc(samples.gamma))$z)       
}else
{
    summary.hyper[3, 1:3] <- c(gamma, gamma, gamma)
    summary.hyper[3, 4:7] <- rep(NA, 4)
}   


summary.results <- rbind(summary.beta, summary.hyper)
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)

    
## Compile and return the results
modelfit <- c(DIC, p.d, WAIC, p.w, LMPL)
names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL")

if(fix.rho.S & fix.rho.T)
{
    samples.rhoext <- NA
}else if(fix.rho.S & !fix.rho.T)
{
    samples.rhoext <- samples.gamma
    names(samples.rhoext) <- "rho.T"
}else if(!fix.rho.S & fix.rho.T)
{
    samples.rhoext <- samples.rho  
    names(samples.rhoext) <- "rho.S"
}else
{
    samples.rhoext <- cbind(samples.rho, samples.gamma)
    colnames(samples.rhoext) <- c("rho.S", "rho.T")
}

if(n.miss==0) samples.Y = NA


samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi),  rho=mcmc(samples.rhoext), tau2=mcmc(samples.tau2), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
model.string <- c("Likelihood model - binomial (logit link function)", "\nLatent structure model - Autoregressive CAR model\n")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL, formula=formula, model=model.string,  X=X)
class(results) <- "carbayesST"
    if(verbose)
    {
    b<-proc.time()
    cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
    return(results)
}
