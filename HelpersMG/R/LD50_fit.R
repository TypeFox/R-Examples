.LD50_fit <- function(par, fixed.parameters=NULL, alive, N, doses, equation) {

  par <- c(par, fixed.parameters)
  p <- .modelLD50(par, doses, equation)
  p <- ifelse(p==0, 1E-9, p)
  p <- ifelse(p==1, 1-1E-9, p)

    if (any(is.infinite(p))) {return(Inf)} else {
   return(-sum(dbinom(alive, N, p, log = TRUE)))
  }
  
}

.modelLD50 <- function(par, doses, equation="logistic") {
#  if (class(par)=="data.frame") par <- na.omit(t(par))[,1]
#  names(parx) <- colnames(par)
#  par <- parx
  # embryogrowth:::.modelTSD(par, doses, equation)
  if (equation=="logistic")	p <- 1/(1+exp((1/par["S"])*(par["P"]-doses)))
  if (equation=="logit")	p <- 1/(1+exp(par["P"]+doses*par["S"]))
  if (equation=="probit")	p <- pnorm(par["P"]+doses*par["S"])
  if (equation=="richards") p <- ifelse(par["K"]>3 & sign(par["P"]-doses)==sign(par["S"]), 
                                        0.5*exp((doses-par["P"])/(par["S"]*exp(par["K"]))), 
                                        (1+(2^exp(par["K"])-1)*exp((1/par["S"])*(par["P"]-doses)))^(-1/exp(par["K"])))
  if (equation=="double-richards") p <- ifelse(doses<par["P"], 
                                               ifelse(par["K1"]>3 & sign(par["P"]-doses)==sign(par["S"]), 
                                                      0.5*exp((doses-par["P"])/(par["S"]*exp(par["K1"]))), 
                                                      (1+(2^exp(par["K1"])-1)*exp((1/par["S"])*(par["P"]-doses)))^(-1/exp(par["K1"]))) ,
                                               ifelse(par["K2"]>3 & sign(par["P"]-doses)==sign(par["S"]), 
                                                      0.5*exp((doses-par["P"])/(par["S"]*exp(par["K2"]))), 
                                                      (1+(2^exp(par["K2"])-1)*exp((1/par["S"])*(par["P"]-doses)))^(-1/exp(par["K2"])))
  )
  if (equation=="hill") if (par["P"]<=0) {p <- rep(Inf, length(doses))} else {p <- 1/(1+exp((1/par["S"])*(log(par["P"])-log(doses))))}
  if (equation=="hulin") {
    elin <- par["K1"]*doses+par["K2"]
    p <- ifelse(elin>6, 0, (1+(2^exp(elin)-1)*exp((1/par["S"])*(par["P"]-doses)))^(-1/exp(elin)))}
  return(p)
  
}
