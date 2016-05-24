GDI <- function (itemBank, item, x, it.given, model=NULL, lower = -4, upper = 4, nqp = 33, 
                type = "GDI", priorDist="norm", priorPar = c(0, 1), D=1, X=NULL, lik = NULL) 
{
if (type != "GDI" & type != "GDIP") 
  stop("'type' must be either 'GDI' or 'GDIP'", call. = FALSE)
if (!is.null(X) & !is.null(lik)) {
  if (length(X) != length(lik)) stop("'X' and 'lik' must have the same length!",call.=FALSE)
}
probs <- NULL
GDIF_1 <-NULL
GDIF_2 <-NULL
par <- rbind(itemBank[item,])
if (is.null(X))
  X <- seq(from=lower,to=upper,length=nqp)
if (is.null(model)) {
  if (is.null(lik)) {
    L <- function(th, r, param) 
          prod(Pi(th, param,D=D)$Pi^r * (1 - Pi(th,param,D=D)$Pi)^(1 - r))
    lik <- sapply(X,L,x,it.given)
  }
  if (type=="GDIP") {
    pd<-switch(priorDist,norm=dnorm(X,priorPar[1],priorPar[2]),unif=dunif(X,priorPar[1],priorPar[2]))
    lik <- lik * pd
  }
  lik <- lik/sum(lik)
  probs[1:nqp] <- Pi(X[1:nqp],par,D=D)$Pi
  GDIF_1[1:nqp] <- probs[1:nqp]^2 
  GDIF_1 <- GDIF_1 * lik
  GDIF_1 <- sum(GDIF_1)
  GDIF_2[1:nqp] <- probs[1:nqp]
  GDIF_2 <- GDIF_2 * lik
  GDIF_2 <- sum(GDIF_2)^2
  crit.value <- GDIF_1 - GDIF_2
}
RES <- crit.value
return(RES)
}
