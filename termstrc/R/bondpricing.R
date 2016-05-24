
## Bond pricing function

bond_prices <- function(method="ns", beta, m, cf, lambda){
  # calculate spot rates
  spot_rates <- spotrates(method,beta,m,lambda)/100

  # replace NaNs by zeros
  spot_rates[is.nan(spot_rates)] <- 0        
     
  # calculate discount factors
  discount_factors <- exp(-m*spot_rates)

  # calculate bond prices
  bond_prices <- apply(cf*discount_factors, 2, sum)  
  
  # return spot rates, discount factors and bond prices
  return(list(spot_rates=spot_rates,discount_factors=discount_factors,
            bond_prices=bond_prices))
}  


## Calculation of bond yields

bond_yields <- function(cashflows, m, searchint=c(-1,1), tol=1e-10) {
  # convert input data to matrices if necessary
  if (!is.matrix(cashflows))
    cashflows <- as.matrix(cashflows)
  if (!is.matrix(m))
    m <- as.matrix(m)
  # create empty bond yields matrix in appropriate size
  bondyields <- matrix(0, nrow=ncol(cashflows), ncol=2)
  # put maturity of every bond into first column of bond yields matrix
  bondyields[,1] <- apply(m, 2, max)
  # traverse list of bonds
  for (i in seq_len(ncol(cashflows))) {
    # present value of cash flows for root finding 
    pvcashflows <- function(y) {
      t(cashflows[,i])%*%exp(-m[,i]*y)
    }
    # calculate bond yields
    bondyields[i,2] <- uniroot(pvcashflows, searchint, tol = tol,maxiter=3000)$root 
  }
  # return calculated bond yields matrix
  rownames(bondyields) <- colnames(cashflows)
  colnames(bondyields) <- c("Maturity","Yield")
  bondyields
}

## Calculation of the duration,modified duration and duration based weights 
duration <- function (cf_p,m_p,y){
  y <- matrix(rep(y,nrow(m_p)),ncol=ncol(m_p),byrow=TRUE)
  # mac cauly duration
  d <- apply(cf_p*m_p*exp(-y*m_p),2,sum)/-cf_p[1,]
  # modified duration
  md <- d/(1+y[1,])
  omega <- (1/d)/sum(1/d)
  dur <- cbind(d,md,omega)
  colnames(dur) <- c("Duration","Modified duration","Weights")
  dur
}

## Loss function: mean squared (weighted error)
loss_function <- function(p,phat,omega) {
  sum(omega*((p-phat)^2))
}
