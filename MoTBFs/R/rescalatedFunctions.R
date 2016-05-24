#' Rescales an MoTBF Function
#' 
#' A collation of function to reescale an MoTBF function
#' to the original limits and scale.
#' 
#' @name rescaledFunctions
#' @rdname rescaledFunctions
#' @param fx A function of class \code{"motbf"} learned from a scaled data.
#' @param data A \code{"numeric"} vector containing the original data without being scaled.
#' @param parameters A \code{"numeric"} vector with the coefficients to create the rescaled MoTBF.
#' @param num A \code{"numeric"} value which contains the denominator of the coefficient
#' in the exponential. By default it is 5.
#' @seealso \link{univMoTBF}
#' @return An \code{"motbf"} function of the original data.
#' @examples
#' ## 1. EXAMPLE
#' X <- rchisq(1000, df = 8) ## data
#' modX <- scale(X) ## scale data
#' 
#' ## Learning
#' f <- univMoTBF(modX, POTENTIAL_TYPE = "MOP", nparam=10) 
#' plot(f, xlim = range(modX), col=2)
#' hist(modX, prob = TRUE, add = TRUE)
#' 
#' ## Rescale
#' origF <- rescaledMoTBFs(f, X) ## rescaledMOP(f, X)
#' plot(origF, xlim = range(X), col=2)
#' hist(X, prob = TRUE, add = TRUE)
#' meanMOP(origF) ## mean(X)
#' 
#' ## 2. EXAMPLE 
#' X <- rweibull(1000, shape = 20, scale= 10) ## data
#' modX <- as.numeric(scale(X)) ## scale data
#' 
#' ## Learning
#' f <- univMoTBF(modX, POTENTIAL_TYPE = "MTE", nparam = 9) 
#' plot(f, xlim = range(modX), col=2, main="")
#' hist(modX, prob = TRUE, add = TRUE)
#' 
#' ## Rescale
#' origF <- rescaledMoTBFs(f, X) ## rescaledMTE(f, X)
#' plot(origF, xlim = range(X), col=2)
#' hist(X, prob = TRUE, add = TRUE)
#' 

#' @rdname rescaledFunctions
#' @export
rescaledMoTBFs <- function(fx, data)
{
  if(is.mop(fx)) f <- rescaledMOP(fx, data)
  if(is.mte(fx)) f <- rescaledMTE(fx, data)
  return(f)
}

#' @rdname rescaledFunctions
#' @export
rescaledMOP <- function(fx, data)
{
  parameters <- coef(fx)
  if(length(parameters)==1) f <- noquote(paste(1/diff(range(data)), "+0*x", sep=""))
  else f <- ToStringRe_MOP(parameters, data) 
  f <- list(Function = f, Subclass = "mop",
            Domain = (fx$Domain)*sd(data)+mean(data))
  f <- motbf(f)
  return(f)
}

#' @rdname rescaledFunctions
#' @export
ToStringRe_MOP <- function(parameters, data) 
{
  str <- parameters[1]/sd(data)
  sign <- parameters; sign[sign<0]=""; sign[sign>=0]="+"
  if(mean(data)<0) smean <- "+" else smean <- "-"
  for(i in 2: length(parameters)) str <- noquote(paste(str,sign[i], parameters[i]/(sd(data)^i), "*(x",smean,abs(mean(data)),")^",i-1, sep=""))
  return(str)  
}

#' @rdname rescaledFunctions
#' @export
rescaledMTE <- function(fx, data)
{
  parameters <- coef(fx)
  parExp <- coeffExp(fx)[-1]
  if(length(parameters)==1) f <- noquote(paste(1/diff(range(data)), "+0*exp(x)", sep=""))
  else f <- ToStringRe_MTE(parameters, data, 1/parExp[1]) 
  f <- list(Function = f, Subclass = "mte",
            Domain = (fx$Domain)*sd(data)+mean(data))
  f <- motbf(f)
  return(f)
}

#' @rdname rescaledFunctions
#' @export
ToStringRe_MTE <- function(parameters, data, num = 5)
{
  mu <- mean(data); sde=sd(data)
  str <- parameters[1]/sde
  sign <- parameters; sign[sign<0]=""; sign[sign>=0] <- "+" 
  if(mu<0) smean <- "+" else smean <- "-"
  
  if((length(parameters)-1)>0) {
    for(i in 2: length(parameters)){
      if(i<4){
        if(i%%2==1) p=-1 else p=1
        if(smean=="-") m <- -mu else m <- mu
        str <- paste(str,sign[i], (parameters[i]/sde)*exp(p*m/(num*sde)), "*exp(",ifelse(i%%2==1, "-", ""),1/(num*sde),"*x",")", sep="")
      }else{
        if(i%%2==1) p <- -(i%/%2) else p <- i%/%2
        if(smean=="-") m <- -mu else m <- mu
        str <- paste(str,sign[i], (parameters[i]/sde)*exp(p*m/(num*sde)), "*exp(",p/(num*sde),"*x)", sep="")
      }    
    }
  }
  else {
    str <- paste(str,"+0*exp(1/num*x)", sep="")
  }
  return(noquote(str))
}

#' @rdname rescaledFunctions
#' @export
meanMOP=function(fx)
{
  fx <- noquote(as.character(fx))
  t <- strsplit(fx, split="*(", fixed = T, perl = FALSE, useBytes = FALSE)[[1]][2]
  t1 <- strsplit(t, split=")", fixed = T, perl = FALSE, useBytes = FALSE)[[1]][1]
  
  mu <- strsplit(t1, split="x", fixed = T, perl = FALSE, useBytes = FALSE)[[1]][2]  
  mu <- as.numeric(mu)
  if(mu<0) mu <- abs(mu) else mu <- as.numeric(paste("-",mu, sep=""))
  return(mu)
}
