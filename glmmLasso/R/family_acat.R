acat <- function()
{
family <- "acat"

  responseFun <- function(eta){
    q <- length(eta)
    eta.help <- matrix(rep(c(0,eta),each=q+1),ncol=q+1)
    eta.help[upper.tri(eta.help)] <- 0
    pi <- cumprod(c(1,exp(eta[-q])))/sum(apply(exp(eta.help),1,prod))
    pi
  } 
    
  linkinv <- function(eta_cat){
    t(apply(eta_cat,1,acat()$responseFun))
  }
  
  createSigmaInv <- function(mu){
    Sigma <- diag(mu) - mu %*% t(mu)
    solve(Sigma)
  }
  
  mulist <- function(mu){
    split(mu, rep(1:nrow(mu), ncol(mu)))
    }
  
  SigmaInv <- function(mu){
    SigmaInv <- as.matrix(bdiag(lapply(acat()$mulist(mu),acat()$createSigmaInv)))
    SigmaInv
  }
  
  createD <- function(mu){
    q <- length(mu)
    
    D2 <- matrix(0,q,q)
    diag(D2) <- -(1/mu)
    
    if(q==2){
      D2[2,1] <- 1/mu[-1]
    }else{
      diag(D2[2:q,1:(q-1)]) <- 1/mu[-1]
    }
    
    D2[,q] <- -1/(1-sum(mu))
    D2[q,q] <- -(1-sum(mu[-q]))/((1-sum(mu))*mu[q])
    
    D <- solve(D2)
    D
  }
  
  deriv.mat <- function(mu_cat){
    as.matrix(bdiag(lapply(acat()$mulist(mu_cat),acat()$createD)))
  }
  
  multivariate <- TRUE
  
  ret.list <- list(responseFun = responseFun, linkinv = linkinv, SigmaInv = SigmaInv,
                   createSigmaInv = createSigmaInv, createD = createD, deriv.mat = deriv.mat,
                   mulist = mulist,
                   multivariate = multivariate, family = family)
  return(ret.list)

}
  

