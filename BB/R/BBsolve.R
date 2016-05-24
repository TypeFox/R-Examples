BBsolve <- function(par, fn, method=c(2,3,1), control=list(), quiet=FALSE, ...) 
    {
    ctrl <- list(maxit = 1500, M = c(50, 10), tol = 1e-07, trace = FALSE, 
        triter = 10, noimp = 100, NM = c(TRUE, FALSE))
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    if(is.matrix(par)) stop("argument par should not be a matrix in BBsolve.")
    ctrl[namc] <- control
    M <- ctrl$M
    maxit <- ctrl$maxit
    tol <- ctrl$tol
    trace <- ctrl$trace
    triter <- ctrl$triter
    noimp <- ctrl$noimp
    NM <- if (length(par) > 1 & length(par) <= 20) ctrl$NM else FALSE

    control.pars <- expand.grid(method=method, M=M, NM=NM)
	
    feval <- iter <-  0
    ans.best.value <- Inf
    for (i in 1: nrow(control.pars) ) {
      cpars <- unlist(control.pars[i, ])
      #cat("Try : ", i, "Method = ", cpars[1], "M = ", cpars[2], "Nelder-Mead = ", cpars[3], "\n")

      temp <- try(dfsane(par=par, fn, method=cpars[1], 
                         control=list(M=as.numeric(cpars[2]), NM=cpars[3], 
			   maxit=maxit, tol=tol, trace=trace, triter=triter, 
			   noimp=noimp),
			 quiet=quiet, alertConvergence=FALSE, ...), silent=TRUE)
      if (!inherits(temp, "try-error")) {
         feval <- feval + temp$feval
         iter <- iter + temp$iter

         if (temp$convergence  == 0) {
              ans.best <- temp
              ans.best$feval <- feval
              ans.best$iter <- iter
              ans.best$cpar <- cpars
              break
              } 
         else if (temp$residual < ans.best.value) {
              ans.best <- temp
              ans.best.value <- ans.best$residual
              ans.best$feval <- feval
              ans.best$iter <- iter
              ans.best$cpar <- cpars
              }
         }
      }  # "i" loop completed

    if(!quiet) {if (ans.best$convergence != 0)
                     cat ("  Unsuccessful convergence.\n")
                else cat ("  Successful convergence.\n")
		}

    ans.best
    }
