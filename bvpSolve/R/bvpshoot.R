## EXTENDED FOR DAEs
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using the shooting method
##==============================================================================

bvpshoot<- function(yini = NULL, x, func, yend = NULL, parms = NULL, 
    order = NULL,  guess = NULL, 
    jacfunc = NULL, bound = NULL, jacbound = NULL, 
    leftbc = NULL, posbound = NULL, ncomp = NULL, 
    atol = 1e-8, rtol = 1e-8, extra = NULL, maxiter = 100, 
    positive = FALSE, method = "lsoda", ...)  {

## ---------------------
## check input
## ---------------------
  
  if (is.null(yini)   && is.null(bound))
    stop("either 'yini' and 'yend' or 'bound' should be inputted")
  if (!is.null(yini)  && is.null(yend))
    stop("if 'yini' is given, 'yend' should also be given")
  if (!is.null(bound) && is.null(posbound)&& is.null(leftbc))
    stop("if 'bound' is given, 'posbound' or 'leftbc' should also be given")
  if (!is.null(yini)  && !is.null(bound))
    stop("either 'yini' or bound should be given, not both")
  if (!is.null(bound)  && !is.null(extra))
    stop("cannot combine 'bound' with 'extra'; use 'yini' and 'yend' ")
  if (is.character(func) | is.character(jacfunc) |
      is.character(bound) | is.character(jacbound))
    stop("cannot use 'bvpshoot' with compiled code")

  lex      <- length(extra)

  # dots for passing to function and dots for method
  fdots <- list(...)
  metargs <- names(formals(method))
  fmet <- which(names(fdots) %in% metargs)
  fdots [fmet] <- NULL

  if (is.function(method)) 
    Ode <- method  else 
    Ode <- function(...) ode (method = method, ...)

  if (length(fdots) == 0)
    Do.call <- function(x, ll) do.call(x, ll)
  else
    Do.call <- function(x, ll) do.call(x, c(ll,unlist(fdots)))

## ---------------------------
## The order of the equations
## ---------------------------
  JacFunc <- jacfunc # save jacfunc for further passing it to ode()
  Func <- func
  testit <- FALSE
  attrib <- rep(0,4) # number of steps, number of fn evaluations, # jacobians ,nr ivp

  if (! is.null(order)) {
    mstar <- sum(order)
    neq   <- length(order)
    if (! is.null(ncomp)) 
      if (sum(order) != ncomp)
        stop("'ncomp' and 'order' not compatible: ncomp should equal sum(order)")
    if (max(order) > 1) { # wrapper over func and jacfunc
      testit <- TRUE
      stareq <- cumsum(order)            # from func to vector returned to c
      higord <- (1:mstar)[-stareq]  # from state to vector returned to c
      
      # expand func
      Fret <- numeric(length = mstar)
      Func    <- function(x, state, parms,...)  {
        FF <- Do.call(func, list(x, state, parms))
        Fret[stareq] <- unlist(FF[1])
        Fret[higord] <- state[higord+1]
        FF[1] <- list(Fret)
        FF
      }

      if (! is.null(jacfunc))
       stop ("can not combine analytical jacobian with higher-order equations - remove 'jacfunc'")
   }
   if (is.null(ncomp)) ncomp <- mstar    
  }             

## ---------------------
## yini or bound
## ---------------------
  if (! is.null(yini))  {    
    # initial value yini; a function or vector
    inity <- function(X, Parms) {  # initialises yini and parms..
      if (is.function(yini))
        Y <- do.call(yini,list(X,Parms))
      else Y <- yini
      if (lini>0)
        Y[inix] <- X[1:lini]
      names(Y)<-Ynames
      Y
    }
    # parameters; some may be estimated (lex)
    initparms <- function (X) {
      Initparms <- parms
      if(lex>0)
        Initparms[1:lex] <- X[(lini+1):(lini+lex)]
      return(Initparms)  
    }  

    if (is.function(yini))
      y <- Do.call(yini,list(extra, parms))
    else
      y <- yini
      
    # root function to solve for  - note jacfunc = NULL to avoid error in devel R 2.12
    rootfun <- function(X, jacfunc = NULL, ...)  {  
      times <- c(x[1], x[length(x)])
      Parms <- initparms(X)
      Y     <- inity(X,Parms)
      out   <- Ode(y=Y, times=times, func=Func, jacfunc=JacFunc,
                 parms=Parms, #method=method,
                 atol=atol, rtol=rtol, ...)
      attrib <<- attrib + c(attributes(out)$istate[c(2,3,14)],1)
    # deviation with yend should be =0             
      if (is.function(yend) )
        Res   <- yend(out[nrow(out),2:(ly+1)], Y, Parms,...)
      else {
        Res <- yend - out[nrow(out),2:(ly+1)]
        Res <- Res[! is.na(Res)]
      }
      return(Res)
    }
    # the jacobian of root function: not specified 
    JacBound <- NULL   

  } else {                  # bound is specified   
    posspecified <- FALSE   
    if (is.null(leftbc)& is.null(posbound))
       stop("'leftbc' or 'posbound' should be inputted if 'bound' is given")
    if (! is.null(posbound)) {
      if (all(posbound %in% c(x[1],x[length(x)])))  
        leftbc <- sum(posbound == x[1])
      else  {
         posspecified <- TRUE 
         # check consistency  
         if (! is.null(ncomp)) {
          if(length(posbound) != ncomp)
            stop("'posbound' should have a length = number of variables ")
         } else ncomp <- length(posbound)
            
         iipos <- which (x %in% posbound)
         if (length(iipos) != ncomp)
           stop("all elements in 'posbound' should also be in 'x' ")
         
      }
    }   
    if (! is.null(guess) & ! is.null(ncomp))
      if (length(guess) != ncomp) stop ("length of 'guess' should be = number of variables")
      
    y <- guess  
  
    rootfun <- function(X, jacfunc = NULL, ...)  {  
      if (! posspecified) 
         times <- c(x[1], x[length(x)])
      else
         times <- x   
      out   <- Ode(y = X, times = times, func = Func, jacfunc = JacFunc,
                   parms = parms, #method = method,
                   atol = atol, rtol = rtol, ...)
      attrib <<- attrib + c(attributes(out)$istate[c(2,3,14)],1)
      Res <- vector(length=ly)
      if (! posspecified) {
        Yend <- out[nrow(out),2:(ly+1)]             
        for (i in 1:leftbc) 
          Res[i] <- Do.call(bound,list(i,X,parms))
        if (leftbc < ly) 
          for (i in  (leftbc+1):ly) 
            Res[i]<- Do.call(bound,list(i,Yend,parms))
      } else 
      for (i in 1:length(posbound)) {
        ii <- iipos[i]          
        Res[i] <- Do.call(bound,list(i,out[ii,-1],parms))
      }
      return(Res)
    }
    JacBound <- NULL   
    if (! is.null(jacbound)) {  
      JAC <- matrix(nrow=length(y),ncol=length(y),0)
      JacBound <- function(x,...)  {
        for (i in 1:ly)    
          JAC[i,]<- Do.call(jacbound,list(i,x,parms))
        return(JAC)  
      }    
    }    
  }
## ---------------------
## names, guess of y
## ---------------------
  
  Ynames <- attr(y,"names")
  if (is.null(Ynames) & ! is.null(yend)) 
    Ynames <- names(yend)
  if (is.null(y)) {
    if (is.null(ncomp ))
      stop("don't know number of variables - provide 'ncomp'")
    y <- rep(0,ncomp)
    warning("estimates for initial conditions not given ('guess'); assuming 0's")
  }  
  if (is.null(y)) 
    stop ("either provide 'guess' for initial conditions or 'ncomp', number of compartments")

  ly <-  length(y)
  
  inix     <- which (is.na(y))
  lini     <- length(inix)
  if (lini > 0 & is.null(guess))  {
    warning("estimates for unknown initial conditions not given ('guess'); assuming 0's")
    guess <- rep(0,lini)
  }

  if (! is.null(yini) & lini != length(guess))  {
    if (is.null(extra))
      stop("length of guess should be equal to number of NAs in y") else
    if (lex > length(parms))
      stop("length of extra should be smaller than number of parameters")
  }
  if (lini > 0)
    y[inix] <- guess

  if (is.null(yini) & is.null(guess)) 
    guess <- y
  if (! is.null(yini) & lini+lex==0)
    stop ("this is not a boundary value problem - use initial value problem solver instead")

  if (testit) {
    if (! is.null(yini))
      Parms <- initparms(extra)
    else
      Parms <- parms
    FF <- Do.call(func,list(x[1], y, Parms))
    if (length(FF[[1]]) != neq)
      stop("function 'func' should return as many elements as the length of 'order'")  
  }
  
## ---------------------
## root solver: 
## ---------------------
  # find unknown initial condition + extra parameter  
  sol <- multiroot(start = c(guess,extra), f = rootfun, atol=atol, rtol=rtol,
                   maxiter=maxiter, jacfunc = JacBound, positive =positive, ...)
  
## ---------------------
## Output 
## ---------------------

  if (! is.null(yini) ) {
    Parms <- initparms(sol$root)
    Y     <- inity(sol$root, Parms)
  }  else {
    Parms <- parms 
    Y     <-  sol$root
  }
    
  out <- Ode (times=x, func=Func, y=Y, parms=Parms, jacfunc=jacfunc, #method=method, 
              atol=atol, rtol=rtol, ...)
  attrib <- attrib + c(attributes(out)$istate[c(2,3,14)],1)

  attrib[2] <- attrib[2] +attrib[4] +1 # correct for one more per ivp solution
              
  attr(out,"roots")  <- data.frame(root=sol$root,
                                   f.root=sol$f.root, iter=sol$iter)
  class(out) <- c("bvpSolve","matrix","deSolve")  # a boundary value problem
  colnames(out)[1] <- "x"
  attr(out,"name") <- "bvpshoot"
  attr(out,"istate2") <- attrib[c(2,3,1,4)]
  names (attr(out,"istate2")) <- c("nfunc", "njac", "nstep", "nivp")

  out
}
