cumulative <- function()
{
family <- "cumulative"

    
  linkinv <- function(eta_cat){
    t(apply(eta_cat,1,binomial()$linkinv))
  }
  
  createSigmaInv <- function(mu){
    Sigma <- mu%*%t(1-mu)
    Sigma[lower.tri(Sigma)]<-t(Sigma)[lower.tri(Sigma)]
    SigmaInv <- try(solve(Sigma),silent=T); k=1
    while(class(SigmaInv)=="try-error")
    {  
    Sigma.new <- Sigma + diag(10^(-7+k),length(mu))
    SigmaInv <- try(solve(Sigma.new),silent=T);k<-k+1
    }
    return(SigmaInv)
  }
  
  mulist <- function(mu){
    split(mu, rep(1:nrow(mu), ncol(mu)))
    }
  
  SigmaInv <- function(mu){
    SigmaInv <- as.matrix(bdiag(lapply(cumulative()$mulist(mu),cumulative()$createSigmaInv)))
    SigmaInv #<- as.matrix(SigmaInv)
  }
  
  deriv.mat <- function(mu_cat){
    eta <- c(apply(mu_cat,1,binomial()$linkfun))
    diag(binomial()$mu.eta(eta))
  }
  
  multivariate <- TRUE
  
  ret.list <- list(linkinv = linkinv, SigmaInv = SigmaInv,
                   createSigmaInv = createSigmaInv, deriv.mat = deriv.mat,
                   mulist = mulist, multivariate = multivariate, family = family)
  return(ret.list)

}
  

