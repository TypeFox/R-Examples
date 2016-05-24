bentlerParameters <-
function(x, N, nFactors, log=TRUE, cor=TRUE,
         minPar=c(min(lambda) - abs(min(lambda)) +.001, 0.001),
         maxPar=c(max(lambda), lm(lambda ~ I(length(lambda):1))$coef[2]),
         resParx=c(0.01, 2), resPary=c(0.01, 2),
         graphic=TRUE, resolution=30, typePlot="wireframe", ...){
 stopMessage  <- paste("\n These indices are only valid with a principal component solution.\n",
                       " ...................... So, only positive eugenvalues are permitted.\n",
                       sep="")
 lambda       <- eigenComputes(x, cor=cor, ...)
 if (length(which(lambda <0 )) > 0) {cat(stopMessage);stop()}
 
 k     <- nFactors
 p     <- length(lambda)
 q     <- p-k
 i     <- 1:q
 x     <- q-i
 l     <- lambda[k+i]
 n     <- N - 1
 
 # Bentler (1996, p. 133) maximization of equations 8 and 9
 f1    <- function(n,l,x,alpha,beta) sum((n*l-(n+1)*(alpha+beta*x))/((alpha+beta*x)^2))
 f2    <- function(n,l,x,alpha,beta) sum((n*l-(n+1)*(alpha+beta*x))*x/((alpha+beta*x)^2))
 f     <- function(alpha,beta) f1(n,l,x,alpha,beta)^2+f2(n,l,x,alpha,beta)^2
 if (log == FALSE)  F <- function(y) f(y[1],y[2])  else  F <- function(y) log(f(y[1],y[2]))
 
 figure <- NULL
 if (graphic == TRUE) {
  p1        <- seq(resParx[1], resParx[2], length=resolution)
  p2        <- seq(resPary[1], resPary[2], length=resolution)
  data      <- expand.grid(Alpha = p1, Beta = p2)
  data      <- data.frame(data, y=numeric(length(data$Alpha)))
  for( i in 1:length(data$Alpha)) data$y[i] <- F(c(data$Alpha[i],data$Beta[i]))

  if (log == FALSE) zlab <- "y" else zlab <- "log(y)"
  if (typePlot == "wireframe")   figure    <- wireframe(  y ~ Alpha * Beta, data=data, zlab=zlab, ...)
  if (typePlot == "contourplot") figure    <- contourplot(y ~ Alpha * Beta, data=data, region=TRUE, ...)
  if (typePlot == "levelplot")   figure    <- levelplot(  y ~ Alpha * Beta, data=data, region=TRUE, ...)
  }
  
 res   <- nlminb(objective=F,start=lm(l~x)$coefficients,lower=c(minPar[1],minPar[2]),upper=c(maxPar[1],maxPar[2]))
 para  <- res$par[1]
 parb  <- res$par[2]
 # Bentler (1996, p. 133) equation 7
 # !!! Warning: Bentler and Yuan (1998) were in error for the definition of LRT !!!
 # !!! So N and n must be inversed in the first logarithm                       !!!
 lrt   <- N*(k-p)*(log(n/N)+1)-N*sum(log(lambda[(k+1):p]/(para+parb*x))) + n*sum(lambda[(k+1):p]/(para+parb*x))
 df    <- q-2
 resp  <- list(convergence=res$convergence, figure=figure, coefficients=res$par,
              lrt=lrt, df=df,k=k,p.value=1-pchisq(lrt,df))
 names(resp$coefficients)<-c("alpha","beta")
 return(resp)
 }