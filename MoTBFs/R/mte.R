#' Fitting Exponential Models
#' 
#' These functions fit mixtures of truncated exponentials (MTEs). 
#' Least square optimization is used to  
#' minimize the quadratic error between the empirical 
#' cumulative distribution and the estimated one. 
#'  
#' @name mte.learning
#' @rdname mte.learning
#' @param X A \code{"numeric"} data vector.
#' @param nparam Number of parameters of the function.
#' @param domain A \code{"numeric"} containing the range where defining the function.
#' @param maxParam A \code{"numeric"} value which indicate the maximum number of coefficients in the function. By default it is \code{NULL}; 
#' if not, the output is the function which gets the best BIC with at most this number of parameters.
#' @return 
#' \code{mte.lerning()} returns a list of n elements:
#' \item{Function}{An \code{"motbf"} object of the \code{'mte'} subclass.}
#' \item{Subclass}{\code{'mte'}.}
#' \item{Domain}{The range where the function is defined to be a legal density function.}
#' \item{Iterations}{The number of iterations that the optimization problem needs to minimize
#' the errors.}
#' \item{Time}{The time which spend the CPU for solving the problem.}
#' 
#' \code{bestMTE()} returns a list including the polynomial function with the best BIC score, 
#' the number of parameters, the best BIC value and an array contained 
#' the BIC values of the evaluated functions.
#' @details 
#' \code{mte.learning()}:
#' The returned value \code{$Function} is the only visible element which contains the mathematical expression. 
#' Using \link{attributes} the name of the others elements are shown and also they can be abstract with \code{$}.
#' The \link{summary} of the function also shows all this elements.
#'
#' \code{bestMTE()}:
#' The first returned value \code{$bestPx} contains the output of the \code{mte.learning()} function
#' with the number of parameters which gets the best BIC value, taking into account the  
#' Bayesian information criterion (BIC) to penalize the functions. It evaluates the two next functions,
#' if the BIC doesn't improve then the function with the last best BIC is returned.
#' 
#' @seealso \link{univMoTBF} A complete function for learning MOPs which includes extra options.
#' @examples
#' 
#' ## 1. EXAMPLE
#' data <- rchisq(1000, df=3)
#' 
#' ## MTE with fix number of parameters
#' fx <- mte.learning(data, nparam=7, domain=range(data))
#' hist(data, prob=TRUE, main="")
#' plot(fx, col=2, xlim=range(data), add=TRUE)
#' 
#' ## Best MTE in terms of BIC
#' fMTE <- bestMTE(data, domain=range(data))
#' attributes(fMTE)
#' fMTE$bestPx
#' hist(data, prob=TRUE, main="")
#' plot(fMTE$bestPx, col=2, xlim=range(data), add=TRUE)
#' 
#' ## 2. EXAMPLE
#' data <- rexp(1000, rate=1/3)
#' 
#'  ## MTE with fix number of parameters
#' fx <- mte.learning(data, nparam=8, domain=range(data))
#' ## Message: The nearest function with odd number of coefficients 
#' hist(data, prob=TRUE, main="")
#' plot(fx, col=2, xlim=range(data), add=TRUE)
#' 
#' ## Best MTE in terms of BIC
#' fMTE <- bestMTE(data, domain=range(data), maxParam=10)
#' attributes(fMTE)
#' fMTE$bestPx
#' attributes(fMTE$bestPx)
#' hist(data, prob=TRUE, main="")
#' plot(fMTE$bestPx, col=2, xlim=range(data), add=TRUE)

#' @rdname mte.learning
#' @export
mte.learning <- function(X, nparam, domain)
{
  ## Time
  tm <- Sys.time()
  ## CDF quadratic function to be minimized
  x <- sort(as.numeric(X)); f <- ecdf(x); y <- f(x)
  y <- unique(y); x <- unique(x);  n <- length(x)
  
  diffDomain <- round(abs(diff(domain))/2)
  seq <- c(0.5,5,50)
  num <- seq[which(abs(seq-diffDomain)==min(abs(seq-diffDomain)))]
  ## N.records is equal 1
  if(length(x)==1){
    if(nparam!=1) return(NULL)
    P <- asMTEString(1/diff(domain), num)
    P <- list(Function = P, Subclass = "mte", Domain = domain,
              Iterations = 0, Time = 0)
    P <- motbf(P)
    return(P)
  }
  
  ## Fit a constant
  if(nparam==1){
    P <- asMTEString(1/diff(domain), num)
    P <- list(Function = P, Subclass = "mte", Domain = domain,
              Iterations = 0, Time = 0)
    P <- motbf(P)
    return(P)
  }
  
  ## Free parameters
  nparam <- nparam-1
  if((nparam%%2)!=0) cat("The nearest function with odd number of coefficients \n")
  if(nparam==1){
    P <- asMTEString(1/diff(domain), num)
    P <- list(Function = P, Subclass = "mte", Domain = domain,
              Iterations = 0, Time = 0)
    P <- motbf(P)
    return(P)
  }
  nparam <- (nparam - nparam%%2)/2
  
  exx <- c(); xx <- c()
  for(i in 1:nparam){
    ex <- exp((i/num)*x)
    nex <- exp((-i/num)*x)
    exx <- rbind(ex, nex)
    xx <- rbind(xx,exx)
  }
  xx <- rbind(x,xx)
  
  ## Dmat
  XX <-  xx%*%t(xx) 
  
  ## dvec
  Xy <-  xx%*%y
  
  ## Compute the constraints
  xnew <- seq(min(domain), max(domain), 10^(-3))
  if(xnew[length(xnew)]!=max(domain)) xnew <- c(xnew,max(domain))
  tt1 <- c()
  ter <- c(); t <- 1
  for(i in 1:nparam){
    ter1 <- (i/num)*exp((i/num)*xnew)
    ter2 <- (-i/num)*exp((-i/num)*xnew)
    ter <- rbind(ter1,ter2)
    t <- rbind(t,ter)
  }
  tt1 <- matrix(c(tt1,t), nrow=2*nparam+1)
  
  tt2 <- c()
  for(i in 1:nparam){
    ex <- exp((i/num)*min(domain))
    nex <- exp((-i/num)*min(domain))
    exx <- rbind(ex, nex)
    tt2 <- rbind(tt2,exx)
  }
  tt2 <- rbind(min(domain),tt2)
  
  tt3 <- c()
  for(i in 1:nparam){
    ex <- exp((i/num)*max(domain))
    nex <- exp((-i/num)*max(domain))
    exx <- rbind(ex, nex)
    tt3 <- rbind(tt3,exx)
  }
  tt3 <- rbind(max(domain),tt3)
  
  ## Amat
  AA <- cbind(tt2, tt3, tt1) 
  
  ## bvec
  B <-  c(0,1, rep(10^(-3), length(xnew)))
  
  ## Solve the optimization problem
  tr <- tryCatch(solve.QP(XX, Xy, AA, B, meq=2), error = function(e) NULL) 
  finaltm <- Sys.time() - tm
  if(is.null(tr)==T){
    return(NULL)
  }else{
    soluc <- tr 
    parameters <- soluc$solution; #parameters
    Px <- asMTEString(parameters, num)
    Px <- list(Function=Px, Subclass="mte")
    Px <- motbf(Px)
    
    ## Derivative (PDF)
    P <- derivMoTBF(Px)
    
    P <- asMTEString(coef(P), num)
    P <- list(Function = P, Subclass = "mte", Domain = domain,
              Iterations = tr$iterations[1], Time = finaltm)
    P <- motbf(P)
    return(P)
  }
}

#' @rdname mte.learning
#' @export
bestMTE <- function(X, domain, maxParam=NULL)
{
  bestBIC <- -10^10; bestfx <- 0; bestter <- 0; bestparam <- 0
  nparam <- 1; i <- 0; MTE1 <- 0; vecBIC <- 0
  if(!is.null(maxParam)){
    Pxs <- list()
    repeat{
      if(!is.null(maxParam)&&(nparam>maxParam)) break
      if(i==4) break
      
      ## Compute an MTE function with a fix number of parameters
      Pp <- mte.learning(X, nparam, domain)
      
      if((is.null(Pp)==T)&&(!is.motbf(bestfx))) {
        nparam <- nparam+2
        i <- i+1
        next
      }
      if((is.null(Pp)==T)&&(is.motbf(bestfx))) break
      
      nparam <- nparam+2
      ## Compute the BIC score
      BiC <- sum(log(as.function(Pp)(X)))
      #BiC <- BICMoTBF(Pp, X)
      vecBIC <- c(vecBIC,BiC)
      Pxs[[length(Pxs)+1]] <- Pp
    }
    p <- which(vecBIC[-1]==max(vecBIC[-1]))
    bestBIC <- vecBIC[p+1]
    bestfx <- Pxs[[p]]
    bestter <- length(coef(bestfx))-2
  }else{
    repeat{
      if(!is.null(maxParam)&&(nparam>maxParam)) break
      if(i==4) break
      
      ## Compute an MTE function with a fix number of parameters
      Pp <- mte.learning(X, nparam, domain)
      
      if((is.null(Pp)==T)&&(!is.motbf(bestfx))) {
        nparam <- nparam+2
        i <- i+1
        next
      }
      if((is.null(Pp)==T)&&(is.motbf(bestfx))) break
      
      ## Compute the BIC score
      BiC <- BICMoTBF(Pp, X)
      vecBIC <- c(vecBIC,BiC); vecBIC
      
      if(length(vecBIC)<=2){
        if(BiC>bestBIC){
          bestBIC <- BiC
          bestfx <- Pp
          bestter <- nparam-2
          nparam <- nparam+2
        } else {
          nparam <- nparam+2
          next
        }
      } else{
        if(BiC>bestBIC){
          bestBIC <- BiC
          bestfx <- Pp
          bestter <- nparam-2
          nparam <- nparam+2
          
        } else{
          if(bestBIC==vecBIC[length(vecBIC)-1]){
            if(length(unique(vecBIC[-1]))==1) break
            nparam <- nparam+2
            next
          } else{
            break
          }
        }
      }
    }
  }
  
  
  result <- list(bestPx=bestfx, bestBIC=bestBIC, paramN=bestter+2, vecBIC=vecBIC[2:length(vecBIC)])
  return(result)
}


#' Parameters to MTE String
#' 
#' This function builds a string with the structure of a \code{'mte'} function.
#' 
#' @param parameters A \code{"numeric"} vector containing the coefficients.
#' @param num A \code{"numeric"} value which contains the denominator of the coefficient
#' in the exponential.
#' @return A \code{"character"} string with an \code{'mte'} structure.
#' @export
#' @examples
#' 
#' param <- -5.8
#' asMTEString(param)
#' 
#' param <- c(5.2,0.3,-3,4)
#' asMTEString(param)
#'  
asMTEString  <-  function(parameters, num = 5) 
{
  str  <-  parameters[1]
  if(length(parameters)==1) return(noquote(paste(str,"+0*exp(", 1/num, "*x)", sep="")))
  sign  <-  parameters; sign[sign<0]=""; sign[sign>=0]="+" 
  for(i in 2: length(parameters)){
    if(i<4) str  <-  paste(str,sign[i], parameters[i], "*exp(",ifelse(i%%2==1, "-", ""),1/num,"*x)", sep="")
    else    str  <-  paste(str,sign[i], parameters[i], "*exp(",ifelse(i%%2==1, "-", ""),(i%/%2)/num,"*x)", sep="")
  }
  return(noquote(str))
}


#' Extract MTE Coefficients
#' 
#' It extracts the parameters of the learned mixtures of truncated exponential models.
#' 
#' @name coef.mte
#' @rdname coef.mte
#' @param fx An \code{"motbf"} function of subclass \code{'mop'}.
#' @return An array with the parameters of the function.
#' @details
#' \code{coeffMOP()} return the coefficients of the terms in the function.
#' 
#' \code{coeffPol()} returns the coefficients of the potential of the polynomial basis in the function. 
#' @seealso \link{coef.motbf} and \link{univMoTBF}
#' @examples
#' 
#' ## 1. EXAMPLE
#' data <- rnorm(1000, mean=5)
#' fx1 <- univMoTBF(data, POTENTIAL_TYPE = "MTE")
#' hist(data, prob=TRUE, main="")
#' plot(fx1, xlim=range(data), col="red", add=TRUE)
#' coeffMTE(fx1) ## coef(fx1)
#' coeffExp(fx1)
#' 
#' ## 2. EXAMPLE
#' data <- rexp(1000, rate=1/2)
#' fx2 <- univMoTBF(data, POTENTIAL_TYPE = "MTE")
#' hist(data, prob=TRUE, main="")
#' plot(fx2, xlim=range(data), col="red", add=TRUE)
#' coeffMTE(fx2) ## coef(fx2)
#' coeffExp(fx2)
#'  

#' @rdname coef.mte
#' @export
coeffMTE <- function(fx)
{
  fx <- noquote(as.character(fx))
  f1 <- substr(fx, 1, 1)
  t <- strsplit(fx, split="-", fixed = T)[[1]]
  for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")##le volvemos a añadir el simbolo negativo
  if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
  t2 <- c()
  for(i in 1:length(t)){
    t1 <- strsplit(t[i], split="+", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
    t2 <- c(t2,t1)
  }
  t3 <- c()
  for(i in 1:length(t2)){
    t1 <- strsplit(t2[i], split="*", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
    t3 <- c(t3,t1)
  }
  pos1 <- grep("e", t3)
  pos <- grep("exp", t3)
  pos1 <- pos1[pos1%in%pos==F]
  if(length(pos1)!=0){
    h1 <- c()
    for(i in pos1){
      if(as.numeric(t3[i+1])>=0) sign="+" else sign=""
      h <- paste(t3[i], sign, t3[i+1], sep="")
      h1 <- c(h1,h)
    }
    
    t3[pos1] <- h1
    t4 <- t3[-(pos1+1)]
    t3 <- t4
  }
  pos <- grep("exp", t3)
  if(t3[1]!="-") parameters <- t3[1] else parameters <- t3[2]
  if(length(pos)!=0) t3 <- t3[pos-1]
  parameters <- c(parameters, t3)
  
  return(as.numeric(parameters))
}

#' @rdname coef.mte
#' @export
coeffExp <- function(fx){
  
  param <- coef(fx)
  fx <- noquote(as.character(fx))
  t <- fx; t2 <- c()
  for(i in 1:length(t)){
    t1 <- strsplit(t[i], split="(", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
    t2 <- c(t2,t1)
  }
  pos <- grep("*x)", t2)
  t3 <- t2[pos]; t2 <- c()
  for(i in 1:length(t3)){
    t1 <- strsplit(t3[i], split="*x)", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
    t2 <- c(t2,t1[1])
  }
  options(warn=-1)
  if(length(param)==(length(t2)+1)) t2 <- c(0,t2)
  return(as.numeric(t2))
}

#' Derivative MTE
#' 
#' Compute the derivative of an \code{"motbf"} object with \code{'mte'} subclass.
#' 
#' @param fx An \code{"motbf"} object of the \code{'mte'} subclass.
#' @return The derivative which is also an \code{"motbf"} function.
#' @seealso \link{univMoTBF} for learning and \link{derivMoTBF} for 
#' general \code{"motbf"} models.
#' @export
#' @examples
#' 
#' ## 1. EXAMPLE
#' X <- rexp(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' derivMTE(Px)
#' 
#' ## 2. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' derivMTE(Px)
#' 
#' ## 3. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' derivMTE(Px)
#' ## Message: It is an 'motbf' function but not 'mte' subclass.
#' class(Px)
#' subclass(Px)
#' 
derivMTE <- function(fx)
{
  if(!is.motbf(fx)) return(cat("It is not an 'motbf' function."))
  if(is.motbf(fx)&&!is.mte(fx)) return(cat("It is an 'motbf' function but not 'mte' subclass."))

  parameters <- coef(fx)
  parExp <- coeffExp(fx)[-1]
  coeffderiv <- c(parameters[1], parameters[-1]*parExp)
  if((length(coeffderiv)==2)&&(coeffderiv[2]==0)) coeffderiv <- rep(0,2)
  P <- asMTEString(coeffderiv, 1/parExp[1])
  P <- list(Function=P, Subclass="mte")
  P <- motbf(P)
  return(P)
}

#' Integral MTE
#' 
#' Method to calculate the non-defined integral of an \code{"motbf"} object of \code{'mte'} subclass.
#' 
#' @param fx An \code{"motbf"} object of subclass \code{'mte'}.
#' @return The non-defined integral of the function.
#' @seealso \link{univMoTBF} for learning and \link{integralMoTBF} 
#' for a more complete function to get defined and non-defined integrals
#' of class \code{"motbf"}.
#' @export
#' @examples
#' 
#' ## 1. EXAMPLE
#' X <- rexp(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' integralMTE(Px)
#' 
#' ## 2. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MTE")
#' integralMTE(Px)
#' 
#' ## 3. EXAMPLE
#' X <- rnorm(1000)
#' Px <- univMoTBF(X, POTENTIAL_TYPE="MOP")
#' integralMTE(Px)
#' ## Message: It is an 'motbf' function but not 'mte' subclass.
#' class(Px)
#' subclass(Px)
#' 
integralMTE <- function(fx)
{  
  if(!is.motbf(fx)) return(cat("It is not an 'motbf' function."))
  if(is.motbf(fx)&&!is.mte(fx)) return(cat("It is an 'motbf' function but not 'mte' subclass"))
  
  parameters <- coeffMTE(fx)
  coefExponential <- coeffExp(fx)[-1]
  str <- paste(parameters[1], "*x", sep="")
  if((length(parameters)-1)>0) {
    for(i in 2:length(parameters)){
      if((parameters[i]*(1/coefExponential[i-1]))>=0) sign <- "+" else sign <- ""
      str <- paste(str,sign, parameters[i]*(1/coefExponential[i-1]), "*exp(",coefExponential[i-1],"*x)", sep="")
    }
  }else {
    ## An MTE constant
    str  <-  paste(str,"+0*exp(", 1/coefExponential[1], "*x)", sep="")
  }
  f <- noquote(str)
  f <- list(Function = f, Subclass = "mte")
  f <- motbf(f)
  return(f)
}

