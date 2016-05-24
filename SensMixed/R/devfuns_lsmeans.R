

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
      
    
      vcov_out <- sigma2 * tcrossprod(envff2$pp$RXi()) 
      
      return(as.matrix(vcov_out))      
    }
  } 
  
  class(ans) <- "vcovJSStheta2.temp"
  ans
}

#update model
updateModel <- function(model, mf.final, reml.lmerTest.private, 
                        l.lmerTest.private.contrast, 
                        devFunOnly.lmerTest.private = FALSE)
{
  if(!mf.final == as.formula(paste(".~.")))
  {
    inds <-  names(l.lmerTest.private.contrast) %in% attr(terms(as.formula(mf.final)), 
                                                          "term.labels")
    #update contrast l.lmerTest.private.contrast
    l.lmerTest.private.contrast <- l.lmerTest.private.contrast[inds]
  }
  
  nfit <- update(object=model, formula.=mf.final, REML=reml.lmerTest.private ,
                 contrasts=l.lmerTest.private.contrast, 
                 devFunOnly = devFunOnly.lmerTest.private, evaluate=FALSE)
  env <- environment(formula(model))
  assign("l.lmerTest.private.contrast", l.lmerTest.private.contrast, envir=env)
  assign("reml.lmerTest.private", reml.lmerTest.private, envir=env)
  assign("devFunOnly.lmerTest.private", devFunOnly.lmerTest.private, envir=env)
  nfit <- eval(nfit, envir = env) 
  return(nfit)   
}

getREML <- function(model)
{
  if(inherits(model,"merMod"))
    return(getME(model, "is_REML"))
  
}



