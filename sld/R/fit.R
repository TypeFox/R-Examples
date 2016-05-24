fit.sld.lmom.given <- function(lmoms,n=NULL){
  if (lmoms[3]>(1/3)) {stop("No QB skew logistic distribution corresponds to these L Moment values.\nThese L Moments are more right skew than the exponential distribution, the limiting case of the QB Skew Logistic.")} 
  if (lmoms[3]<(-1/3)) {stop("No QB skew logistic distribution corresponds to these L Moment values.\nThese L Moments are more left skew than the reflected exponential distribution, the limiting case of the QB Skew Logistic.")}
  ah = lmoms[1] - 6*(lmoms[3]*lmoms[2]) # alpha hat 
  bh = 2*lmoms[2] # beta hat
  dh = 0.5*(1+3*lmoms[3]) # delta hat
  lmomest <- c(ah,bh,dh)
  names(lmomest) <- c("alpha","beta","delta")
  if (!is.null(n)){ # Sample size is known - calculate Std Errors
    om = dh*(1-dh) # omega
    se.alpha = bh * sqrt((57 + (125*pi^2-1308)*om)/(15*n))
    se.beta = bh * sqrt(4/(3*n) * (1 - (pi^2-8)*om))
    se.delta = sqrt((8-(397+160*om-20*pi^2*(om+2))*om)/(15*n))
    lmomse <- c(se.alpha,se.beta,se.delta)
    ret <- cbind(lmomest,lmomse)
    dimnames(ret) <- list(c("alpha","beta","delta"),
        c("Estimate","Std. Error"))
  } else {
    ret <- lmomest # return just the estimates
  }
  ret
}

fit.sld.lmom <- function(data){
  fit.sld.lmom.given(lmom.sample(data,3),n=length(data))
}

lmom.sample <- function(data,max.mom=3){
  LM <- samlmu(x=data,nmom=max.mom)
  LM
}