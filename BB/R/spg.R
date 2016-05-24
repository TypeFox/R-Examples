spg <- function(par, fn, gr=NULL, method=3, lower=-Inf, upper=Inf, 
	   project=NULL, projectArgs=NULL, 
	   control=list(), quiet=FALSE, alertConvergence=TRUE,  ... ) {


  box <- if (any(is.finite(upper)))   TRUE
    else if (any(is.finite(lower)))   TRUE
    else                              FALSE  


  prj <- if (box)                     TRUE
    else if (!is.null(project))       TRUE
    else                              FALSE  

  if (is.character(project)) project <- get(project, mode="function")
 
  if (box){
    if (is.null(project) ){
         # upper and lower for default. Expand if scalar
         if(is.null(projectArgs)) projectArgs <- list()
	 
	 if( (!is.null(projectArgs$lower)) | (!is.null(projectArgs$upper))) 
	    warning("Using lower and upper spg arguments, ", 
	            "not using those specified in projectArgs.")

         projectArgs$lower <- 
	     if (length(lower)==1) rep(lower, length(par)) else lower
         projectArgs$upper <- 
	     if (length(upper)==1) rep(upper, length(par)) else upper
         
         # default previously called projectBox
	 project <-  function(par, lower, upper) {
             # Projecting to ensure that box-constraints are satisfied
             par[par < lower] <- lower[par < lower]
             par[par > upper] <- upper[par > upper]
             return(par)
             }
         }

    if (identical(project, projectLinear)){
	 
       if( (!is.null(projectArgs$lower)) | (!is.null(projectArgs$upper))) 
	    warning("Using lower and upper spg arguments, ", 
	            "not using those specified in projectArgs.")

       if(is.null(projectArgs$A)) stop(
	  "projectLinear requires the A matrix to be specified in projectArgs.")

       if(is.null(projectArgs$b)) stop(
	  "projectLinear requires the b vector to be specified in projectArgs.")

        # upper and lower. Expand if scalar
       if (length(lower)==1) lower <- rep(lower, length(par))
       if (any(zi <- is.finite(lower))){
	  projectArgs$A <- rbind(projectArgs$A, diag(length(par))[zi,])
	  projectArgs$b <-     c(projectArgs$b, lower[zi])
	  }
 
       if (length(upper)==1) upper <- rep(upper, length(par))
       if (any(zi <- is.finite(upper))){
	  projectArgs$A <- rbind(projectArgs$A, diag(-1, length(par))[zi,])
	  projectArgs$b <-     c(projectArgs$b, -upper[zi])
	  }
       }
  }


  # control defaults
  # Added `ftol' to the control list:  RV change on 02-06-2011 
  ctrl <- list(M=10, maxit=1500, ftol=1.e-10, gtol=1.e-05, maxfeval=10000, maximize=FALSE, 
        trace=TRUE, triter=10, quiet=FALSE, eps=1e-7, checkGrad=NULL, checkGrad.tol=1.e-06) 
  namc <- names(control)
  if (! all(namc %in% names(ctrl)) )
     stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])     

  ctrl[namc ] <- control
  M	   <- ctrl$M
  maxit    <- ctrl$maxit
  gtol     <- ctrl$gtol
  ftol     <- ctrl$ftol     # RV change on 02-06-2011
  maxfeval <- ctrl$maxfeval
  maximize <- ctrl$maximize
  trace    <- ctrl$trace
  triter   <- ctrl$triter
  eps      <- ctrl$eps
  checkGrad <- ctrl$checkGrad
  checkGrad.tol <- ctrl$checkGrad.tol  
    
  grNULL <- is.null(gr)  
  fargs <- list(...)


  func <- if (maximize) function(par, ...) c(-fn(par, ...))
                   else function(par, ...) c( fn(par, ...))
  
  # first evaluate the function to be sure the initial guess works
  # and use the timing to decide if the analytic gradient should be checked.
  # c() in next is for case of a 1x1 matrix value

  f.time <- system.time(f <- try(func(par, ...),silent=TRUE))
  feval <- 1
  
  # set the default for checking the gradient based on how long
  # a function evaluation takes, and the number of parameters.
  # 6* is a very crude approximation of how many funtion evaluations are
  # necessary for each dimension.
  if (is.null(checkGrad))  
     if  (((f.time[1]+f.time[2]) * 6*length(par)) < 10) {
        checkGrad <-TRUE 
        }
     else {
        checkGrad <- FALSE
        if(!grNULL) warning(
	 "Default checking of gradient turned off because of time require.",
	 "See the help for spg to enable this.")
	} 

  if (class(f)=="try-error" )
        stop("Failure in initial function evaluation!", f)
  else if ( !is.numeric(f) || 1 != length(f) )
        stop("function must return a scalar numeric value!")
  else if (is.nan(f) | is.infinite(f) | is.na(f) )
        stop("Failure in initial function evaluation!")
    
  f0 <- fbest <- f
 

  
  ################ local function
  nmls <- function(p, f, d, gtd, lastfv, feval, func, maxfeval, fargs ){
    # Non-monotone line search of Grippo with safe-guarded quadratic interpolation
    gamma <- 1.e-04
    fmax <- max(lastfv)
    alpha <- 1
    pnew <- p + alpha*d
    fnew <- try(do.call(func, append(list(pnew) , fargs )),silent=TRUE)
    feval <- feval + 1
 
    if (class(fnew)=="try-error" | is.nan(fnew) )
            return(list(p=NA, f=NA, feval=NA, lsflag=1))
 
    while(fnew > fmax + gamma*alpha*gtd) {
    	if (alpha <= 0.1) alpha <- alpha/2
    	else {
    	    atemp <- -(gtd*alpha^2) / (2*(fnew - f - alpha*gtd))
    	    if (atemp < 0.1 | atemp > 0.9*alpha) atemp <- alpha/2
    	    alpha <- atemp
    	    }

    	pnew <- p + alpha*d
    	fnew <- try(do.call(func, append(list(pnew), fargs )), silent=TRUE)
    	feval <- feval + 1
 
    	if (class(fnew)=="try-error" | is.nan(fnew) )
	       return(list(p=NA, f=NA, feval=NA, lsflag=1))
    	if (feval > maxfeval)
	       return(list(p=NA, f=NA, feval=NA, lsflag=2))
 
    	}  #while condition loop ends
 
    return(list(p=pnew, f=fnew, feval=feval, lsflag=0))
    }
  #############################################
  if (!grNULL & checkGrad) {
    requireNamespace("numDeriv", quietly = TRUE)
    grad.num <- numDeriv::grad(x=par, func=fn, ...) 
    grad.analytic <- gr(par, ...)
    max.diff <- max(abs((grad.analytic - grad.num) / (1 + abs(fn(par, ...)))))
    if(!max.diff < checkGrad.tol) {
      cat("Gradient check details:  max. relative difference in gradients= ",
	         max.diff,
	   "\n\n  analytic gradient:",  grad.analytic,
	   "\n\n  numerical gradient:", grad.num
	   )
      stop("Analytic gradient does not seem correct! See comparison above. ",
           "Fix it, remove it, or increase checkGrad.tol." )
      }
    }
  ################ local function 
  # Simple gr numerical approximation. Using func, f and eps from calling env.  	
  # used when user does not specify gr.
  if (grNULL) gr  <-function(par, ...) {
    	df <- rep(NA,length(par))
    	for (i in 1:length(par)) {
    	  dx <- par
    	  dx[i] <- dx[i] + eps 
    	  df[i] <- (func(dx, ...) - f)/eps
    	 }
    	df
	}

  #############################################

  #  Initialization
  lmin <- 1.e-30
  lmax <- 1.e30
  iter <-  0
  lastfv <- rep(-1.e99, M)
  fbest <- NA
  fchg <- Inf          # RV change on 02-06-2011
 
  # this switch is not needed for the numerical grad because the
  #  sign of func is switched.
  grad <- if (maximize & !grNULL) function(par, ...) -gr(par, ...)
                    else          function(par, ...)  gr(par, ...)


  # Project initial guess
  if (prj){
     par <- try(do.call(project,  append(list(par), projectArgs)), silent=TRUE)
 
     if (class(par) == "try-error") 
           stop("Failure in projecting initial guess!", par)
     }
  if (any(is.nan(par), is.na(par)) ) stop("Failure in initial guess!")
  pbest <- par
 
  g <- try(grad(par, ...),silent=TRUE)
 
  if (class(g)=="try-error" ) 
        stop("Failure in initial gradient evaluation!", g)
  else if (any(is.nan(g)) ) 
        stop("Failure in initial gradient evaluation!")
 
  lastfv[1] <- fbest <- f
 
  pg <- par - g
  if (prj){
    pg <- try(do.call(project,  append(list(pg), projectArgs)),silent=TRUE)
    if (class(pg)=="try-error" ) stop("Failure in initial projection!", pg)
    } 
    
  if (any(is.nan(pg))) stop("Failure in initial projection!")
 
  pg <- pg - par

  pg2n <- sqrt(sum(pg*pg))
  pginfn <- max(abs(pg))
  gbest <- pg2n
  if (pginfn != 0) lambda <- min(lmax, max(lmin, 1/pginfn))
 
  if (trace) cat("iter: ",0, " f-value: ", f0 * (-1)^maximize, " pgrad: ",pginfn, "\n")

  #######################
  #  Main iterative loop
  #######################
  lsflag <- NULL # for case when tol is already ok initially and while loop is skipped
  while( pginfn > gtol & iter <= maxit & fchg > ftol) {     # RV change on 02-06-2011 
      iter <- iter + 1
 
   d <- par - lambda * g
   if (prj){
     d <- try(do.call(project,  append(list(d), projectArgs)), silent=TRUE)
     if (class(d) == "try-error" | any(is.nan(d))  ) {
          lsflag <- 4
          break
          }
        }
 
      d <- d - par
      gtd <- sum(g * d)
 
      if(is.infinite(gtd)){
        lsflag <- 4
        break
        }
 
      nmls.ans <- nmls(par, f, d, gtd, lastfv, feval , func, maxfeval, fargs)
      lsflag <- nmls.ans$lsflag
 
      if(lsflag != 0) break

      fchg <- abs(f - nmls.ans$f)  # RV change on 02-06-2011 
      f     <- nmls.ans$f
      pnew  <- nmls.ans$p
      feval <- nmls.ans$feval
      lastfv[(iter %% M) + 1] <- f
 
      gnew <- try(grad(pnew, ...),silent=TRUE)     
 
      if (class(gnew)=="try-error" | any(is.nan(gnew)) ){
        lsflag <- 3
        break
        }
 
      s <- pnew - par
      y <- gnew - g
      sts <- sum(s*s)
      yty <- sum(y*y)
      sty <- sum(s*y)
 
      if (method==1) lambda <- 
	   if (sts==0  | sty < 0)  lmax else min(lmax, max(lmin, sts/sty))
      else 
      if (method==2) lambda <- 
           if (sty < 0 | yty == 0) lmax else min(lmax, max(lmin, sty/yty))
      else 
      if (method==3) lambda <-
           if (sts==0  | yty == 0) lmax else min(lmax, max(lmin, sqrt(sts/yty)))
 
 
      par <- pnew
      g   <- gnew
 
      pg <- par - g
      if (prj){
        pg <- try(do.call(project, append(list(pg), projectArgs)), silent=TRUE)
        if (class(pg) == "try-error" | any(is.nan(pg)) ) {
  	  lsflag <- 4
  	  break
  	  }
	}

      pg <- pg - par
      pg2n <- sqrt(sum(pg*pg))
      pginfn <- max(abs(pg))
 
      f.rep <- (-1)^maximize * f
      if (trace && (iter%%triter == 0))
           cat("iter: ",iter, " f-value: ", f.rep, " pgrad: ",pginfn, "\n")
 
      if (f < fbest) {
  	fbest <- f
  	pbest <- pnew
  	gbest <- pginfn
  	}
 
      }   # while condition loop concludes

  if (is.null(lsflag)) {
        if (!quiet) warning("convergence tolerance satisified at intial parameter values.")
	lsflag <- 0
	}
 
  if (lsflag==0) {
    if (pginfn <= gtol | fchg <= ftol) conv <- list(type=0, message="Successful convergence")
    if (iter >= maxit)  conv <- list(type=1, message="Maximum number of iterations exceeded")
    f.rep <- (-1)^maximize * fbest  # This bug was fixed by Ravi Varadhan.  March 29, 2010.
    par <- pbest
    } else {
      par <- pbest
      f.rep <- f <- (-1)^maximize * fbest
      pginfn <- gbest
      if (lsflag==1) conv <- list(type=3, message="Failure:  Error in function evaluation")
      if (lsflag==2) conv <- list(type=2, message="Maximum function evals exceeded")
      if (lsflag==3) conv <- list(type=4, message="Failure:  Error in gradient evaluation")
      if (lsflag==4) conv <- list(type=5, message="Failure:  Error in projection")
      }
 
  if(alertConvergence && ( 0 != conv$type))
          warning("Unsuccessful convergence.")

  return(list(par=par, value=f.rep, gradient =pginfn, 
      fn.reduction=(-1)^maximize * (f0 - f), 
      iter=iter, feval=feval, convergence=conv$type, message=conv$message))
  }
 
 
 
 
