# Carrier prob conditional on the phenotype
carrierprobpheno <- function(fit, method="data", mode="dominant", q=0.02)
{
  
  theta <- fit$parms.est 
  base.dist <- attr(fit, "base.dist")
  agemin <- attr(fit, "agemin")
  data <- attr(fit, "data")
  
  if(sum(is.na(data$mgene))==0) stop("Mutatioin carrier status are all known")
  
  beta.sex <- theta[3]
  beta.gen <- theta[4]
  
  xbeta1 <- beta.sex*data$gender+beta.gen*1
  xbeta0 <- beta.sex*data$gender+beta.gen*0
  
  y <- data$time-agemin
  delta <- data$status
  haz <-hazards(base.dist, y, theta)
  Haz <-cumhaz(base.dist, y, theta)
  
  haz1 <- haz*exp(xbeta1)
  haz0 <- haz*exp(xbeta0)
  Haz1 <- Haz*exp(xbeta1)
  Haz0 <- Haz*exp(xbeta0)

  p <- data$carrp
  if(is.null(p)) p <- carrierprobgeno(data=data, method=method, mode=mode, q=q)$carrp
  
  p1 <- (haz1^delta)*exp(-Haz1)*p
  p0 <- (haz0^delta)*exp(-Haz0)*(1-p)
  
  carrp.pheno <- p1/(p1+p0)
  carrp.pheno[!is.na(data$mgene)] <- data$mgene[!is.na(data$mgene)]
 
  data$carrp.pheno <- carrp.pheno
  
  return(data)

}