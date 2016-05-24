### 'simRasch' : simulate the data matrix from a Rasch Model
###
### Author: Zhushan 'Mandy' Li
### Last update: Sept 06, 2006
### Update Notes: Not finished yet
###       Right Now it can only simulate 1D Rasch Model

## require(mvtnorm)

simRasch <- function(ncat, nitem, nexaminee, beta=NULL){
  simRasch.1D(ncat, nitem, nexaminee, beta)
}

simRasch.1D <- function(ncat, nitem, nexaminee,  beta=NULL)
  #item.mtx=NULL,, trait.cov=NULL)
{
  ## the unidimensional, multicategory,  IRT function
  pIRF <- function(theta=0, a, b){
    wts <- exp(b + a*theta);
    sum <- apply(wts, 1, sum);
    prob <- sweep(wts,1,sum, FUN="/");
    return(prob);
  }


  ## input a, b
  parm <- data.frame(item=1:nitem)
  parm$a <- matrix(rep(0:(ncat-1), each=nitem), ncol=ncat)
  if(is.null(beta)){
    beta <- matrix(rnorm(nitem*(ncat-1), mean=0, sd=1), ncol=(ncat-1))
  }
  parm$b <- cbind(rep(0,nitem), beta)
  
  
  ## Simulation

  ## Draw theta from N(0,1)
  theta <- rnorm(nexaminee);  ## draw theta form N(0,1)
  sim.data <- matrix(NA, nrow=nexaminee, ncol=nitem);

  for(j in (1:nexaminee)) {
    prob <- pIRF(theta[j], parm$a, parm$b);
    cumprob <- t(apply(prob, 1, cumsum))
    jju <- runif(nitem);
    sim.data[j,] <-   apply(jju > cumprob,1, sum)
  }

  result <- list(data=sim.data, beta=beta, theta=theta)
  return(result)
}
