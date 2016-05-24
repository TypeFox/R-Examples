gaussian.CARanova <- function(formula, data=NULL, W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, fix.rho.S=FALSE, rho.S=NULL, fix.rho.T=FALSE, rho.T=NULL, verbose=TRUE)
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
which.miss <- as.numeric(!is.na(Y))
which.miss.mat <- matrix(which.miss, nrow=K, ncol=N, byrow=FALSE)
n.miss <- N.all - sum(which.miss)
Y.miss <- Y
Y.miss[which.miss==0] <- median(Y, na.rm=TRUE)    
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
X.short <- X.standardised[which.miss==1, ]    


#### Offset variable
offset <- try(model.offset(frame), silent=TRUE)
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,N.all)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    
    

        
    
#### Specify the initial parameter values
mod.glm <- glm(Y~X.standardised-1, offset=offset)
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

res.temp <- Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=0, sd = res.sd)
delta <- rnorm(n=N, mean=0, sd = res.sd)
tau2.phi <- var(phi)/10
tau2.delta <- var(delta)/10
nu2 <- runif(1, 0, res.sd)

if(fix.rho.S)
{
    rho <- rho.S
}else
{
    rho <- runif(1)       
}

if(fix.rho.T)
{
    lambda <- rho.T
}else
{
    lambda <- runif(1)       
}   

    
    
#### Check and specify the priors
## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
    if(is.null(prior.nu2)) prior.nu2 <- c(0.001, 0.001)

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
    
    if(length(prior.nu2)!=2) stop("the prior value for nu2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.nu2)) stop("the prior value for nu2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.nu2))!=0) stop("the prior value for nu2 has missing values.", call.=FALSE)       



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
if(lambda<0 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
if(lambda>1 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  


#### Set up matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.delta <- array(NA, c(n.keep, N))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 2))
colnames(samples.tau2) <- c("tau2.phi", "tau2.delta")    
if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
if(!fix.rho.T) samples.lambda <- array(NA, c(n.keep, 1))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.like <- array(NA, c(n.keep, N.all))
samples.deviance <- array(NA, c(n.keep, 1))
if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

  

#### Specify the Metropolis quantities
accept.all <- rep(0,4)
accept <- accept.all
proposal.sd.rho <- 0.02
proposal.sd.lambda <- 0.02
tau2.phi.shape <- prior.tau2[1] + K/2
tau2.delta.shape <- prior.tau2[1] + N/2
nu2.shape <- prior.nu2[1] + N*K/2    
    
    
    
#### .S quantities
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
    
    
    
    #### .T quantities
    ## .T neighbourhood matrix
    D <-array(0, c(N,N))
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(abs((i-j))==1)  D[i,j] <- 1 
        }    
    }
    
    
    ## Create the triplet object
    D.triplet <- c(NA, NA, NA)
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(D[i,j]>0)
            {
                D.triplet <- rbind(D.triplet, c(i,j, D[i,j]))     
            }else{}
        }
    }
    D.triplet <- D.triplet[-1, ]     
    D.n.triplet <- nrow(D.triplet) 
    D.triplet.sum <- tapply(D.triplet[ ,3], D.triplet[ ,1], sum)
    D.neighbours <- tapply(D.triplet[ ,3], D.triplet[ ,1], length)
    
    
    ## Create the start and finish points for W updating
    D.begfin <- array(NA, c(N, 2))     
    temp <- 1
    for(i in 1:N)
    {
        D.begfin[i, ] <- c(temp, (temp + D.neighbours[i]-1))
        temp <- temp + D.neighbours[i]
    }
    
    
    ## Create the determinant     
    if(!fix.rho.T) 
    {
        Dstar <- diag(apply(D,1,sum)) - D
        Dstar.eigen <- eigen(Dstar)
        Dstar.val <- Dstar.eigen$values
        det.Q.D <-  0.5 * sum(log((lambda * Dstar.val + (1-lambda))))    
    }else
    {} 
       
#### Specify quantities that do not change
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE)
Y.mat.trans <- t(Y.mat)
Y.mat.miss <- matrix(Y.miss, nrow=K, ncol=N, byrow=FALSE)
Y.mat.trans.miss <- t(Y.mat.miss)
which.miss.mat.trans <- t(which.miss.mat)
ntime.miss <- apply(which.miss.mat,1,sum)
nspace.miss <- apply(which.miss.mat,2,sum)
phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)
delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)
    


#### Beta update quantities
data.precision.beta <- t(X.short) %*% X.short
    if(length(prior.var.beta)==1)
    {
    prior.precision.beta <- 1 / prior.var.beta
    }else
    {
    prior.precision.beta <- solve(diag(prior.var.beta))
    }



## Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
if(rho==1) tau2.phi.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   
if(lambda==1) tau2.delta.shape <- prior.tau2[1] + 0.5 * (N-1)   



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
        ##################
        ## Sample from nu2
        ##################
        nu2.offset <- as.numeric(Y.mat - offset.mat - regression.mat - phi.mat - delta.mat)[which.miss==1]
        nu2.scale <- prior.nu2[2]  + sum(nu2.offset^2)/2
        nu2 <- 1 / rgamma(1, nu2.shape, scale=(1/nu2.scale)) 
        
        
        
        ####################
        ## Sample from beta
        ####################
        fc.precision <- prior.precision.beta + data.precision.beta / nu2
        fc.var <- solve(fc.precision)
        beta.offset <- as.numeric(Y.mat - offset.mat - phi.mat - delta.mat)[which.miss==1]
        beta.offset2 <- t(X.short) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
        fc.mean <- fc.var %*% beta.offset2
        chol.var <- t(chol(fc.var))
        beta <- fc.mean + chol.var %*% rnorm(p)        
        regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  
        
        
        ####################
        ## Sample from phi
        ####################
        phi.offset <- Y.mat - offset.mat - regression.mat -  delta.mat
        phi.offset2 <- apply(phi.offset,1, sum, na.rm=TRUE)
        temp1 <- gaussiancarupdate(W.triplet, W.begfin, W.triplet.sum, K, phi, tau2.phi, nu2, phi.offset2, rho, ntime.miss)
        phi <- temp1
        if(rho<1)
        {
            phi <- phi - mean(phi)
        }else
        {
            phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
        phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)    
        
        
        
        ####################
        ## Sample from delta
        ####################
        delta.offset <- Y.mat - offset.mat - regression.mat -  phi.mat
        delta.offset2 <- apply(delta.offset,2, sum, na.rm=TRUE)
        temp2 <- gaussiancarupdate(D.triplet, D.begfin, D.triplet.sum, N, delta, tau2.delta, nu2, delta.offset2, lambda, nspace.miss)
        delta <- temp2
        delta <- delta - mean(delta)
        delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)

                
        
        #######################
        ## Sample from tau2.phi
        #######################
        temp2.phi <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, rho)
        tau2.phi.scale <- temp2.phi + prior.tau2[2] 
        tau2.phi <- 1 / rgamma(1, tau2.phi.shape, scale=(1/tau2.phi.scale))
        
        
        #########################
        ## Sample from tau2.delta
        #########################
        temp2.delta <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, lambda)
        tau2.delta.scale <- temp2.delta + prior.tau2[2] 
        tau2.delta <- 1 / rgamma(1, tau2.delta.shape, scale=(1/tau2.delta.scale))
        
        
        
        ##################
        ## Sample from rho
        ##################
        if(!fix.rho.S)
        {
        proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)   
        temp3 <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, proposal.rho)
        det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
        logprob.current <- det.Q.W - temp2.phi / tau2.phi
        logprob.proposal <- det.Q.proposal - temp3 / tau2.phi
        prob <- exp(logprob.proposal - logprob.current)
        
        #### Accept or reject the proposal
        if(prob > runif(1))
        {
            rho <- proposal.rho
            det.Q.W <- det.Q.proposal
            accept[1] <- accept[1] + 1           
        }else
        {
        }              
        accept[2] <- accept[2] + 1           
        }else
        {}
        
        
        #####################
        ## Sample from lambda
        #####################
        if(!fix.rho.T)
        {
        proposal.lambda <- rtrunc(n=1, spec="norm", a=0, b=1, mean=lambda, sd=proposal.sd.lambda)   
        temp3 <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, proposal.lambda)
        det.Q.proposal <- 0.5 * sum(log((proposal.lambda * Dstar.val + (1-proposal.lambda))))              
        logprob.current <- det.Q.D - temp2.delta / tau2.delta
        logprob.proposal <- det.Q.proposal - temp3 / tau2.delta
        prob <- exp(logprob.proposal - logprob.current)
        
        #### Accept or reject the proposal
        if(prob > runif(1))
        {
            lambda <- proposal.lambda
            det.Q.D <- det.Q.proposal
            accept[3] <- accept[3] + 1           
        }else
        {
        }              
        accept[4] <- accept[4] + 1           
        }else
        {}
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        fitted <- as.numeric(offset.mat + regression.mat + phi.mat  + delta.mat)
        deviance.all <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),N.all), log=TRUE)
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
            samples.delta[ele, ] <- delta
            if(!fix.rho.S) samples.rho[ele, ] <- rho
            if(!fix.rho.T) samples.lambda[ele, ] <- lambda
            samples.nu2[ele, ] <- nu2
            samples.deviance[ele, ] <- deviance
            samples.fitted[ele, ] <- fitted
            samples.tau2[ele, ] <- c(tau2.phi, tau2.delta)
            samples.like[ele, ] <- like
            if(n.miss>0) samples.Y[ele, ] <- rnorm(n=n.miss, mean=fitted[which.miss==0], sd=sqrt(nu2))
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
            accept.rho <- 100 * accept[1] / accept[2]
            if(is.na(accept.rho)) accept.rho <- 45
            accept.lambda <- 100 * accept[3] / accept[4]
            if(is.na(accept.lambda)) accept.lambda <- 45
            accept.all <- accept.all + accept
            accept <- rep(0,4)
            
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
            #### lambda tuning parameter
            if(accept.lambda > 50)
            {
                proposal.sd.lambda <- min(proposal.sd.lambda + 0.1 * proposal.sd.lambda, 0.5)
            }else if(accept.lambda < 40)              
            {
                proposal.sd.lambda <- proposal.sd.lambda - 0.1 * proposal.sd.lambda
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
if(!fix.rho.S)
{
    accept.rho <- 100 * accept.all[1] / accept.all[2]
}else
{
    accept.rho <- NA    
}

if(!fix.rho.T)
{
    accept.lambda <- 100 * accept.all[3] / accept.all[4]
}else
{
    accept.lambda <- NA    
}
accept.final <- c(rep(100,3), accept.rho, accept.lambda)
names(accept.final) <- c("beta", "phi", "delta", "rho.S", "rho.T")        
    
   

## Compute DIC
median.phi <- apply(samples.phi, 2, median)
median.delta <- apply(samples.delta, 2, median)  
median.phi.mat <- matrix(rep(median.phi, N), byrow=F, nrow=K)
median.delta.mat <- matrix(rep(median.delta, K), byrow=T, nrow=K)
median.beta <- apply(samples.beta,2,median)
regression.mat <- matrix(X.standardised %*% median.beta, nrow=K, ncol=N, byrow=FALSE)   
fitted.median <- as.numeric(offset.mat + regression.mat + median.phi.mat + median.delta.mat)    
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),N.all), log = TRUE), na.rm=TRUE)
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
    CPO[j] <- 1/median((1 / dnorm(Y[j], mean=samples.fitted[ ,j], sd=sqrt(samples.nu2))))    
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
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

summary.hyper <- array(NA, c(5, 7))     
rownames(summary.hyper) <- c("tau2.S", "tau2.T", "nu2", "rho.S", "rho.T")   
summary.hyper[1,1:3] <- quantile(samples.tau2[ ,1], c(0.5, 0.025, 0.975))
summary.hyper[2,1:3] <- quantile(samples.tau2[ ,2], c(0.5, 0.025, 0.975))
summary.hyper[3,1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,1])), geweke.diag(mcmc(samples.tau2[ ,1]))$z)     
summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,2])), geweke.diag(mcmc(samples.tau2[ ,2]))$z)   
summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.nu2)), geweke.diag(mcmc(samples.nu2))$z)  

if(!fix.rho.S)
{
    summary.hyper[4,1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[4, 4:7] <- c(n.keep, accept.rho, effectiveSize(mcmc(samples.rho)), geweke.diag(mcmc(samples.rho))$z)  
}else
{
    summary.hyper[4, 1:3] <- c(rho, rho, rho)
    summary.hyper[4, 4:7] <- rep(NA, 4)
}

if(!fix.rho.T)
{
    summary.hyper[5, 1:3] <- quantile(samples.lambda, c(0.5, 0.025, 0.975))
    summary.hyper[5, 4:7] <- c(n.keep, accept.lambda, effectiveSize(samples.lambda), geweke.diag(samples.lambda)$z)
}else
{
    summary.hyper[5, 1:3] <- c(lambda, lambda, lambda)
    summary.hyper[5, 4:7] <- rep(NA, 4)
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
    samples.rhoext <- samples.lambda
    names(samples.rhoext) <- "rho.T"
}else if(!fix.rho.S & fix.rho.T)
{
    samples.rhoext <- samples.rho  
    names(samples.rhoext) <- "rho.S"
}else
{
    samples.rhoext <- cbind(samples.rho, samples.lambda)
    colnames(samples.rhoext) <- c("rho.S", "rho.T")
}
if(n.miss==0) samples.Y = NA
colnames(samples.tau2) <- c("tau2.S", "tau2.T")

samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi),  delta=mcmc(samples.delta), tau2=mcmc(samples.tau2), nu2=mcmc(samples.nu2), rho=mcmc(samples.rhoext), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))        
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nLatent structure model - spatial and temporal main effects\n")
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
