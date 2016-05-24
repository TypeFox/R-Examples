poisson.localisedCAR <- function(formula, data=NULL, G, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.delta = NULL, verbose=TRUE)
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
                   
 

#### Response variable
## Create the response
Y <- model.response(frame)

## Check for errors
if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
n <- length(Y)
int.check <- n-sum(ceiling(Y)==floor(Y))
if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)
log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
which.miss <- as.numeric(!is.na(Y))
n.miss <- n - sum(which.miss)

#### Offset variable
## Create the offset
offset <- try(model.offset(frame), silent=TRUE)

## Check for errors
if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
if(is.null(offset))  offset <- rep(0,n)
if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)



#### Design matrix
## Create the matrix
X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
ptemp <- ncol(X)

if(ptemp==1)
{
    X <- NULL
    regression.vec <- rep(0, n)
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
    mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    regression.vec <- X.standardised %*% beta
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
  
     

#### Priors
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

if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    

if(is.null(prior.delta)) prior.delta <- 10
if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    

     
#### Initial parameter values
res.temp <- log.Y - regression.vec - offset
Z <- sample(1:G, size=n, replace=TRUE)
lambda <- sort(runif(G, min=min(res.temp), max=max(res.temp)))
delta <- runif(1,1, min(2, prior.delta))
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=n, mean=rep(0,n), sd = res.sd)
for(i in 1:G)
{
    phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
}
tau2 <- var(phi) / 10

     
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


     
     
## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.phi <- array(NA, c(n.keep, n))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.Z <- array(NA, c(n.keep, n))
samples.lambda <- array(NA, c(n.keep, G))
samples.delta <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))
samples.like <- array(NA, c(n.keep, n))
samples.fitted <- array(NA, c(n.keep, n))

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
proposal.sd.phi <- 0.1
proposal.sd.delta <- 0.1
proposal.sd.lambda <- 0.01
tau2.posterior.shape <- prior.tau2[1] + 0.5 * n



      
     
#### CAR quantities
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(!is.symmetric.matrix(W)) stop("W is not symmetric.", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    


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
     if(!is.null(X))
     {
         proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
         proposal.beta <- beta
         offset.temp <- phi + offset + lambda[Z]  
         for(r in 1:n.beta.block)
         {
             proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
             prob <- poissonbetaupdate(X.standardised, n, p, beta, proposal.beta, offset.temp, Y, prior.mean.beta, prior.var.beta, which.miss)
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
     }else{}
     
                 
               
     ####################
     ## Sample from phi
     ####################
     phi.offset <- regression.vec + offset  + lambda[Z]
     temp1 <- poissoncarupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=n, phi=phi, tau2=tau2, y=Y, phi_tune=proposal.sd.phi, rho=1, offset=phi.offset, which.miss)
     phi <- temp1[[1]]
          for(i in 1:G)
          {
          phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
          }
     accept[1] <- accept[1] + temp1[[2]]
     accept[2] <- accept[2] + n    
         
         
     ##################
     ## Sample from tau2
     ##################
     temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, n, phi, phi, 1)
     tau2.posterior.scale <- temp2 + prior.tau2[2] 
     tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
          
         
    ####################
    ## Sample from lambda
    ####################
    proposal <- c(-1000, lambda, 1000)
          for(i in 1:G)
          {
           proposal[(i+1)] <- rtrunc(n=1, spec="norm", a=proposal[i], b=proposal[(i+2)], mean=proposal[(i+1)], sd=proposal.sd.lambda)    
          }
     proposal <- proposal[2:(G+1)]
     lp.current <- lambda[Z] + phi + regression.vec + offset
     lp.proposal <- proposal[Z] + phi + regression.vec + offset
     prob1 <- sum((exp(lp.current) - exp(lp.proposal)))
     prob2 <- sum(Y * (lp.proposal - lp.current))
     prob <- exp(prob1 + prob2)
          if(prob > runif(1))
          {
          lambda <- proposal
          accept[3] <- accept[3] + 1  
          }else
          {
          }
     accept[4] <- accept[4] + 1       

     
    
    ################
    ## Sample from Z
    ################
    Z.offset <- phi + offset + regression.vec
    Z.proposal <- sample(1:G, size=n, replace=TRUE)
    prior <- delta * ((Z - Gstar)^2 - (Z.proposal-Gstar)^2)    
    like <- exp(Z.offset) * (exp(lambda[Z]) - exp(lambda[Z.proposal])) + Y * (lambda[Z.proposal] - lambda[Z])         
    prob <- exp(like + prior)   
    test <- prob> runif(n)         
    Z[test] <- Z.proposal[test]         

         
         
    ##################
    ## Sample from delta
    ##################
    proposal.delta <-  rtrunc(n=1, spec="norm", a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)    
    prob1 <- sum((Z-Gstar)^2) * (delta - proposal.delta)        
    prob2 <- n * log(sum(exp(-delta *(1:G - Gstar)^2))) - n * log(sum(exp(-proposal.delta *(1:G - Gstar)^2)))
    prob <- exp(prob1 + prob2)    
          if(prob > runif(1))
          {
          delta <- proposal.delta
          accept[5] <- accept[5] + 1  
          }else
          {
          }
     accept[6] <- accept[6] + 1       

         
         
    #########################
    ## Calculate the deviance
    #########################
     lp <- lambda[Z] + phi + regression.vec + offset
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
        samples.phi[ele, ] <- phi
        samples.lambda[ele, ] <- lambda
        samples.tau2[ele, ] <- tau2
        samples.Z[ele, ] <- Z
        samples.delta[ele, ] <- delta
        samples.deviance[ele, ] <- deviance
        samples.like[ele, ] <- like
        samples.fitted[ele, ] <- fitted
        
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
         accept.phi <- 100 * accept[1] / accept[2]
        accept.lambda <- 100 * accept[3] / accept[4]
        accept.delta <- 100 * accept[5] / accept[6]

        
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
         #### lambda tuning parameter
            if(accept.lambda > 40)
            {
            proposal.sd.lambda <- proposal.sd.lambda + 0.1 * proposal.sd.lambda
            }else if(accept.lambda < 20)              
            {
            proposal.sd.lambda <- proposal.sd.lambda - 0.1 * proposal.sd.lambda
            }else
            {
            }              

        #### delta tuning parameter
            if(accept.delta > 50)
            {
            proposal.sd.delta <- min(proposal.sd.delta + 0.1 * proposal.sd.delta, prior.delta/6)
            }else if(accept.delta < 40)              
            {
            proposal.sd.delta <- proposal.sd.delta - 0.1 * proposal.sd.delta
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
accept.phi <- 100 * accept.all[1] / accept.all[2]
accept.lambda <- 100 * accept.all[3] / accept.all[4]
accept.delta <- 100 * accept.all[5] / accept.all[6]
accept.tau2 <- 100

if(!is.null(X))
{
    accept.beta <- 100 * accept.all[7] / accept.all[8]   
    accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.tau2)
    names(accept.final) <- c("beta", "lambda", "delta", "phi", "tau2")   
}else
{
    accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.tau2)
    names(accept.final) <- c("lambda", "delta", "phi", "tau2")   
}



     
## Deviance information criterion (DIC)
median.phi <- apply(samples.phi, 2, median)
median.Z <- round(apply(samples.Z,2,median),0)
median.lambda <- apply(samples.lambda,2,median)
    if(!is.null(X))
    {
    median.beta <- apply(samples.beta,2,median)
    regression.vec <- as.numeric(X.standardised %*% median.beta)   
    }else
    {
    regression.vec <- rep(0,n)
    }

fitted.median <- exp(regression.vec  + median.lambda[median.Z] + median.phi + offset)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.median, log=TRUE))
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
LMPL <- sum(log(CPO))  
  

     
 
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
summary.hyper <- array(NA, c(2 ,7))
summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)
summary.hyper[2, 1:3] <- quantile(samples.delta, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:7] <- c(n.keep, accept.delta, effectiveSize(samples.delta), geweke.diag(samples.delta)$z)

samples.lambda <- mcmc(samples.lambda)
summary.lambda <- t(apply(samples.lambda, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(samples.lambda), geweke.diag(samples.lambda)$z)
Z.used <- as.numeric(names(table(samples.Z)))
summary.lambda <- summary.lambda[Z.used, ]

if(!is.null(X))
{
    samples.beta.orig <- mcmc(samples.beta.orig)
    summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    rownames(summary.beta) <- colnames(X)
    summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)  
    row.names(summary.results)[(p+1):nrow(summary.results)] <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
    colnames(summary.results) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}else
{
    summary.results <- rbind(summary.lambda, summary.hyper)
    row.names(summary.results) <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
    colnames(summary.results) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}



#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, median)
residuals <- as.numeric(Y) - fitted.values

     
## Compile and return the results
modelfit <- c(DIC, p.d, WAIC, p.w, LMPL)
names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL")  
       
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects  model - Localised CAR model\n")
if(is.null(X)) samples.beta.orig = NA

samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi), lambda=mcmc(samples.lambda), Z=mcmc(samples.Z), tau2=mcmc(samples.tau2), delta=mcmc(samples.delta), fitted=mcmc(samples.fitted))          
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=median.Z,  formula=formula, model=model.string, X=X)
class(results) <- "carbayes"

     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
return(results)
}
