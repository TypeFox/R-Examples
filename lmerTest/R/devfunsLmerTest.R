## the modified devfun functions of the lme4 package
## in order to calculate hessian at the optima theta
## the modified vcov function - in order to calculate gradient at the optima theta


## modified devfun3: now depends on theta and not covariance parameters
devfun5 <- function (fm,  reml = TRUE) 
{
  stopifnot(is(fm, "merMod"))
  
  np <- length(fm@pp$theta)
  nf <- length(fixef(fm)) 
  if (!isGLMM(fm)) 
    np <- np + 1L
  n <- nrow(fm@pp$V)
  

  
  ff <- updateModel(fm, .~., getREML(fm), 
                    attr(model.matrix(fm),"contrasts"), 
                    devFunOnly.lmerTest.private = TRUE) 
  
  envff <- environment(ff)
  
  if (isLMM(fm)) {
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      
      #.Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")
      ff(thpars[-np])
      
      sigsq <- thpars[np]^2
      dev <- envff$pp$ldL2() + (envff$resp$wrss() + envff$pp$sqrL(1))/sigsq + n * 
        log(2 * pi * sigsq)      
      if(reml){
        p <- ncol(envff$pp$RX())
        dev <- dev + 2*determinant(envff$pp$RX())$modulus - p * log(2 * pi * sigsq)              
      }
      
      return(dev)     
    }
  }
  
  attr(ans, "thopt") <- fm@pp$theta
  class(ans) <- "devfun5"
  ans
}



## not calling the C code
vcovJSStheta2 <- function(fm)
{
  stopifnot(is(fm, "merMod"))
  
  np <- length(fm@pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  

  
  ff2 <- updateModel(fm, .~., getREML(fm), 
                     attr(model.matrix(fm),"contrasts"), 
                    devFunOnly.lmerTest.private = TRUE) 
  
  envff2 <- environment(ff2)
  
  if (isLMM(fm)) {
    ans <- function(Lc, thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      
      sigma2 <- thpars[np]^2
      ff2(thpars[-np])
      
      
      #.Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(envff2$pp$RXi()) 
      
      return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))        
    }
  } 
  
  class(ans) <- "vcovJSStheta2"
  ans
}


## not calling the C code
vcovJSStheta2.temp <- function(fm)
{
  stopifnot(is(fm, "merMod"))
  
  #vlist <- sapply(fm@cnms, length)
  
  #pp <- fm@pp$copy()
  
  
  #resp <- fm@resp$copy()
  np <- length(fm@pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  
  
  
  ff2 <- updateModel(fm, .~., getREML(fm), 
                     attr(model.matrix(fm),"contrasts"), 
                     devFunOnly.lmerTest.private = TRUE) 
  
  envff2 <- environment(ff2)
  
  if (isLMM(fm)) {
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      
      sigma2 <- thpars[np]^2
      ff2(thpars[-np])
      
      
      #.Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(envff2$pp$RXi()) 
      
      return(as.matrix(vcov_out))      
    }
  } 
  
  class(ans) <- "vcovJSStheta2.temp"
  ans
}

