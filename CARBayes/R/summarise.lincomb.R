summarise.lincomb <- function(model, columns=NULL, quantiles=0.5, distribution=FALSE)
{
#### Summarise a linear combination of the covariates
X <- model$X
n <- nrow(X)
p <- ncol(X)               
beta <- model$samples$beta
n.sample <- nrow(beta)
     
     if(!is.numeric(columns)) stop("the columns argument is not numeric.", call.=FALSE)
N.all <- length(columns)
     if(sum(floor(columns)==ceiling(columns)) < N.all) stop("the columns argument is not all integers.", call.=FALSE)
     if(min(columns)<1) stop("some of the columns chosen are less than 1.", call.=FALSE)
     if(max(columns)>p) stop("some of the columns chosen are greater than the size of the samples object.", call.=FALSE)

     if(!is.numeric(quantiles)) stop("the quantiles argument is not numeric.", call.=FALSE)
     if(min(quantiles) < 0) stop("some of the quantiles are less than zero.", call.=FALSE)
     if(max(quantiles) > 1) stop("some of the quantiles are greater than one.", call.=FALSE)

     if(!is.logical(distribution)) stop("The distribution argument must be true or false.", call.=FALSE)
     
#### Compute the posterior distribution
beta.samples <- beta[ ,columns]
samples <- array(NA, c(n.sample, n))
     for(i in 1:n)
     {
     samples[ ,i] <- beta.samples %*% X[i,columns]     
     }
     
samples.quantiles <- apply(samples, 2, quantile, quantiles) 
    
     
#### Create the results list
     if(distribution==TRUE)
     {
     results <- list(quantiles=t(samples.quantiles), posterior=samples)     
     }else
     {
     results <- list(quantiles=t(samples.quantiles), posterior=NULL)             
     }
}