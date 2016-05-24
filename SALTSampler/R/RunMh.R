RunMh <- function(center, B, concentration = 1, h, type = 'user', 
                  dat = NULL, pars = NULL) {
  #Run Metropolis Hasting constrained to a simplex.
  #User can define a target distribution or can use
  #built-in functions for a multinomial or dirichlet
  
  #Check inputs
  if (!exists("center")) {
    stop("center is not defined")
  }
  if (!is.numeric(center)){
    stop("center is not numeric")
  }
  if (!is.vector(center)){
    stop("center is not a vector")
  }
  
  if (!exists("B")) {
    stop("B is not defined")
  }
  if (!is.vector(B)){
    stop("B is not a vector")
  }
  if (B%%1 != 0){
    stop("B is not an integer")
  }
  if (length(B) != 1){
    stop("B is not of length 1")
  }
  
  if (!exists("concentration")) {
    stop("concentration is not defined")
  }
  if (!is.vector(concentration)){
    stop("concentration is not a vector")
  }
  if (length(concentration) != 1){
    stop("concentration is not of length 1")
  }
  
  if (!exists("h")){
    stop("h is not defined")
  }
  if (!is.numeric(h)){
    stop("h is not numeric")
  }
  if (!is.vector(h)){
    stop("h is not a vector")
  }
  
  if (!exists("type")) {
    stop("type is not defined")
  }
  if (!(type == 'dirichlet' || type == 'user' || type == 'multinom')) {
    stop("Type not recognized. Use 'dirichlet', 'multinom' or 'user'.")
  }
  
  if (!exists("dat")) {
    stop("dat is not defined, although its value can be set to NULL")
  }
  if (!is.null(dat) & (!is.matrix(dat) & !is.vector(dat))){
    stop("dat is not a matrix")
  }
  
  if (!exists("pars")) {
    stop("pars is not defined")
  }
  if (!is.null(pars) & !is.list(pars)){
    stop("pars is not a list")
  }
  
  if (length(center) != length(h)){
    stop("Length of center does not equal length of h")
  }
  
  if((type == 'multinom') & concentration != 1){
    warning("Concentration is specifed for type = 'multinom'. Concentration will be ignored")
  }
  
  if((type == 'multinom') & (is.null(dat))){
    stop("Type 'multinom' selected, but dat is set to NULL")
  }
  
  #Check user-defined Target distribution function (if type = 'user')
  if (type == 'user'){
    
    if (!exists("Target")) {
      stop("Type 'user' is selected, but no target function is defined in the environment")
    }    
    
    if (length(formals("Target")) != 5) {
      stop("Incorrect number of arguments in user-defined Target function.
           Target function must take in the exact argument sequence: ycand, ycurrent, a, dat, pars")
    }
    
    if (any(!names(formals(Target)) == c("ycand", "ycurrent", "a", "dat", "pars"))){
      stop("Incorrect arguments and/or incorrect ordering of arguments in user-defined Target function. 
           Target function must take in the exact argument sequence: ycand, ycurrent, a, dat, pars")
    }
    
    tryCatch (Target(ycand = Logit(c(0.2, 0.3, 0.7)), ycurrent = Logit(c(0.2, 0.3, 0.7)), a = 1,
                     dat = dat, pars = pars), 
              error = function(err) {
                cat("User-defined Target function fails on the following test case:", "\n",
                    "ycand = Logit(c(0.2, 0.3, 0.5)), ycurrent = Logit(c(0.2, 0.3, 0.5)), a = 1,dat = dat, pars = pars)", "\n",
                    "It is producing the following error: ", "\n")
                err$message
              })
    
    testTarget <-  Target(ycand = Logit(c(0.2, 0.3, 0.7)), ycurrent = Logit(c(0.2, 0.3, 0.5)), a = 1,
                          dat = dat, pars = pars)
    if (length(testTarget) != 1) {
      stop("User-defined Target function produces output that does not have length one")
    }
    if (!is.vector(testTarget)) {
      stop("User-defined Target function produces output that is not a vector")
    }
    if (!is.numeric(testTarget)) {
      stop("User-defined Target function produces non-numeric output")
    }
    }
  
  #Select target distribution function
  #(if type = 'multinom' or type = 'dirichlet)
  if (type == 'dirichlet') {
    Target <- function(ycand, ycurrent, a, dat = NULL, pars = NULL) {
      out <- sum((a - 1)*(LogPq(ycand)$logp - LogPq(ycurrent)$logp))
      return(out)
    }
  } else if (type == 'multinom') {
    Target <- function(ycand, ycurrent, a, dat = NULL, pars = NULL) {
      out <- sum(apply(dat, 1, function(x) {x*LogPq(ycand)$logp})) - 
        sum(apply(dat, 1, function(x) {x*LogPq(ycurrent)$logp}))
      return(out)
    }
  }
  
  #Redefining parameters, setting initial time, and making
  #empty vectors and matrices
  zz=proc.time();
  p <- length(center)
  Y <- array(0, c(B, p)) 
  moveCount <- rep(0, p)
  center <- center/sum(center)
  a <- concentration*center 
  ycurrent <- Logit(center)
  
  #Run sampler
  for (i in 1:B) { 
    for (j in sample(p)) { 
      
      #Propose new y
      ycand <- PropStep(ycurrent, j, h[j])
      
      #Decide to accept or reject
      if (any(is.na(ycand) | is.infinite(ycand) | is.nan(ycand))) {
        move <- FALSE 
      } else {
        move <- (log(runif(1)) < attr(ycand, 'dbt') + 
                   Target(ycand, ycurrent, a, dat, pars))
      }
      if (!is.na(move) & move) {
        ycurrent <- ycand
        moveCount[j] <- moveCount[j] + 1
      }
    }
    
    #Store y for this iteration
    Y[i, ] <- ycurrent
  }
  
  #Timing
  runTime <- proc.time()-zz
  
  #Transform samples to natural scale
  S <- matrix(exp(LogPq(Y)$logp), nrow(Y))
  
  #Return results
  return(list(Y = Y, S = S, runTime = runTime, moveCount = moveCount, p = p, 
              center = center, B = B, concentration = concentration, h = h, 
              type = type, dat = NULL, a = a))
  }
