## =============================================================================
## print the return code settings - all except rk and daspk
## =============================================================================

printidid <- function(idid) {
  cat(paste("\n  return code (idid) = ", idid), "\n")

  if (idid == 2 || idid ==0)  cat("  Integration was successful.\n") else
  if (idid == 3)  cat("  Integration was successful and a root was found before reaching the end.\n") else
  if (idid == -1) cat("  Excess work done on this call. (Perhaps wrong Jacobian type MF.)\n") else
  if (idid == -2) cat("  Excess accuracy requested. (Tolerances too small.)\n") else
  if (idid == -3) cat("  Illegal input detected. (See printed message.)\n") else
  if (idid == -4) cat("  Repeated error test failures. (Check all input.)\n") else
  if (idid == -5) cat("  Repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.)\n") else
  if (idid == -6) cat("  Error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)\n") else
  if (idid == -7) cat("  Work space insufficient to finish (see messages).\n") else
  if (idid == -8) cat("  A fatal error came from sparse solver CDRV by way of DPRJS or DSOLSS.\n")
}

## =============================================================================
## print the return code settings - all except rk and daspk
## =============================================================================

printidid_rk <- function(idid) {
  cat(paste("\n  return code (idid) = ", idid), "\n")

  if (idid == 2 || idid ==0)  cat("  Integration was successful.\n") else
  if (idid == -1) cat("  Excess work done on this call. (Perhaps maxstep exceeded.)\n") else
  if (idid == -2) cat("  Excess accuracy requested. (Tolerances too small.)\n") else
    cat("  Unknown error code, please inform package developers.\n")
}



## =============================================================================
## print the return code settings - only daspk
## =============================================================================

printidid_daspk <- function(idid) {
  cat(paste("\n  return code (idid) = ", idid), "\n")
    if (idid > 0)   {
      cat ("  integration was succesful\n")
      if (idid == 1) cat("  A step was successfully taken in the intermediate-output mode.  The code has not yet reached TOUT.\n")
      if (idid == 2) cat("  The integration to TSTOP was successfully completed (T = TSTOP) by stepping exactly to TSTOP.\n")
      if (idid == 3) cat("  The integration to TOUT was successfully completed (T = TOUT) by stepping past TOUT. Y(*) and YPRIME(*) are obtained by interpolation.\n")
      if (idid == 4) cat("  The initial condition calculation, with INFO(11) > 0, was successful, and INFO(14) = 1. No integration steps were taken, and the solution is not considered to have been started.\n")
    } else if (idid < 0 & idid > -33)  {
      cat ("  integration was interrupted\n")
      if (idid == -1) cat("  A large amount of work has been expended (about 500 steps).\n") else
      if (idid == -2) cat("  The error tolerances are too stringent.\n") else
      if (idid == -3) cat("  The local error test cannot be satisfied because a zero component in ATOL was specified and the corresponding computed solution component is zero.  Thus, a pure relative error test is impossible for this component.\n") else
      if (idid == -5) cat("  There were repeated failures in the evaluation or processing of the preconditioner (in jacfunc).\n") else
      if (idid == -6) cat("  DDASPK had repeated error test failures on the last attempted step.\n") else
      if (idid == -7) cat("  The nonlinear system solver in the time integration could not converge.\n") else
      if (idid == -8) cat("  The matrix of partial derivatives appears to be singular (direct method).\n") else
      if (idid == -9) cat("  The nonlinear system solver in the time integration failed to achieve convergence, and there were repeated error test failures in this step.\n") else
      if (idid == -10) cat("  The nonlinear system solver in the time integration failed to achieve convergence because IRES was equal to -1.\n") else
      if (idid == -11) cat("  IRES = -2 was encountered and control is being returned to the calling program.\n") else
      if (idid == -12) cat("  DDASPK failed to compute the initial Y, YPRIME.\n") else
      if (idid == -13) cat("  Unrecoverable error encountered inside user's PSOL routine, and control is being returned to the calling program.\n") else
      if (idid == -14) cat("  The Krylov linear system solver could not achieve convergence.\n")
    } else if (idid ==-33)  {
      cat ("  integration was terminated\n")
      cat("  The code has encountered trouble from which it cannot recover.  A message is printed explaining the trouble and control is returned to the calling program.\n")
    }
}

## =============================================================================
## print the integer diagnostics
## =============================================================================

printIstate <- function(istate, name, all = TRUE) {
df <- c( "The return code :",                                              #1
         "The number of steps taken for the problem so far:",                #2
         "The number of function evaluations for the problem so far:",       #3
         "The number of Jacobian evaluations so far:",                       #4
         "The method order last used (successfully):",                       #5
         "The order of the method to be attempted on the next step:",        #6
         "If return flag =-4,-5: the largest component in error vector",     #7
         "The length of the real work array actually required:",             #8
         "The length of the integer work array actually required:",          #9
         "The number of matrix LU decompositions so far:",                   #10
         "The number of nonlinear (Newton) iterations so far:",              #11
         "The number of convergence failures of the solver so far ",         #12
         "The number of error test failures of the integrator so far:",      #13
         "The number of Jacobian evaluations and LU decompositions so far:", #14,
         "The method indicator for the last succesful step,
           1=adams (nonstiff), 2= bdf (stiff):" ,                            #15
         "The current method indicator to be attempted on the next step,
           1=adams (nonstiff), 2= bdf (stiff):",                             #16
         "The number of nonzero elements in the sparse Jacobian:" ,           #17
         "The order (or maximum order) of the method:",                       #18
         "The number of convergence failures of the linear iteration so far", #19
         "The number of linear (Krylov) iterations so far ",                  #20
         "The number of psol calls so far:")                                  #21
  if (name =="mebdfi")
    df[19:21] <- c(
         "The number of backsolves so far",
         "The number of times a new coefficient matrix has been formed so far",
         "The number of times the order of the method has been changed so far")

#  if (is.na(istate[14])) istate[14]<-istate[4]+istate[10]  # Jacobian+LU
  cat("\n--------------------\n")
  cat("INTEGER values\n")
  cat("--------------------\n")
  if (all) ii <- 1:19 else ii <- which(!is.na(istate))
  printmessage(df[ii], istate[ii], Nr=ii)
}

## =============================================================================
## print the real diagnostics
## =============================================================================

printRstate <- function( rstate) {
  if(is.null(rstate)) return()
  df <- c( "The step size in t last used (successfully):",
    "The step size to be attempted on the next step:",
    "The current value of the independent variable which the solver has reached:",
    "Tolerance scale factor > 1.0 computed when requesting too much accuracy:",
    "The value of t at the time of the last method switch, if any:")

  cat("--------------------\n")
  cat("RSTATE values\n")
  cat("--------------------\n")

  ii <- which(!is.na(rstate))

  printmessage(df[ii], rstate[ii])

}

## =============================================================================
## print all diagnostic messages
## =============================================================================

diagnostics.deSolve <- function(obj, Full = FALSE, ...) {
  Attr <- attributes(obj)
  name <- Attr$type

  istate <- Attr$istate
  rstate <- Attr$rstate

    cat("\n--------------------\n")
    cat(paste(name,"return code"))
    cat("\n--------------------\n")

  idid <- istate[1]
  if (name == "lsodes" && idid == -7) idid <- -8
  if (name == "rk")  printidid_rk
  if (name == "daspk") printidid_daspk(idid) else  printidid(idid)

  printIstate(istate, name, all=Full)

  if (name != "rk") printRstate(rstate)

  if (!is.null(Attr$nroot)) {
    cat("--------------------\n")
    cat("ROOT + event \n")
    cat("--------------------\n")
      cat("\n root found at times :",
        signif(Attr$troot, digits = 5), "\n")
  }

  if (name == "lsodar" ||
      (name %in% c("lsode","lsodes","radau") && !is.null(Attr$iroot))) {
    cat("--------------------\n")
    cat("ROOT\n")
    cat("--------------------\n")
    iroot <- which (Attr$iroot ==1)
    if (length (iroot) > 0)  {
      cat("\n root found for root equation:",
        signif(iroot, digits = 0), "\n")
      cat("\n at time :",
        signif(Attr$troot, digits = 5), "\n")

      }
    else
     if (is.null(Attr$nroot)) cat("\n NO root found \n")
    invisible(list(istate=istate, rstate=rstate, iroot = iroot))

  } else
  invisible(list(istate=istate, rstate=rstate))
}

diagnostics.default <- function(obj, ...)
  warning("No diagnostics available for class '", class(obj), "'")

diagnostics <- function(obj, ...) UseMethod("diagnostics")
