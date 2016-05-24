lbfgsb3 <- function(prm, fn, gr=NULL, lower = -Inf, upper = Inf,
         control=list(), ...){

# ?? need to add controls
# if (is.null(gr)) require(numDeriv) # eventually change to "grfwd" etc.
# interface to Fortran Lbfgsb.3.0
## 150121 There is currently no limit on function or gradient evaluations

    tasklist <- c('NEW_X', 'START', 'STOP', 'FG',  # 1-4
       'ABNORMAL_TERMINATION_IN_LNSRCH', 'CONVERGENCE', #5-6
       'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL',#7
       'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH',#8
       'ERROR: FTOL .LT. ZERO', #9
       'ERROR: GTOL .LT. ZERO' ,#10
       'ERROR: INITIAL G .GE. ZERO', #11
       'ERROR: INVALID NBD', # 12
       'ERROR: N .LE. 0', # 13
       'ERROR: NO FEASIBLE SOLUTION', # 14
       'ERROR: STP .GT. STPMAX', # 15
       'ERROR: STP .LT. STPMIN', # 16
       'ERROR: STPMAX .LT. STPMIN', # 17
       'ERROR: STPMIN .LT. ZERO', # 18
       'ERROR: XTOL .LT. ZERO', # 19
       'FG_LNSRCH', # 20
       'FG_START', # 21
       'RESTART_FROM_LNSRCH', # 22
       'WARNING: ROUNDING ERRORS PREVENT PROGRESS', # 23
       'WARNING: STP .eq. STPMAX', # 24
       'WARNING: STP .eq. STPMIN', # 25
       'WARNING: XTOL TEST SATISFIED')# 26
# CONV in 6, 7, 8; ERROR in 9-19; WARN in 23-26
 
# if (!is.loaded("lbfgsb3.so")) dyn.load("lbfgsb3.so") # get the routines attached

# control defaults -- idea from spg
ctrl <- list(maxit = 500, trace = 0, iprint = 0L)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control


# Here expand control list, but for moment leave alone
      iprint <- as.integer(ctrl$iprint)
      factr <- 1.0e+7
      pgtol <- 1.0e-5
      nmax <- 1024L
      mmax <- 17L
if (length(prm) > nmax) stop("The number of parameters cannot exceed 1024")
      n <- as.integer(length(prm))
      m <- 5L # default 
## Define the storage
nbd <- rep(2L, n) # start by defining them "on" -- adjust below
nwa<-2*mmax*nmax + 5*nmax + 11*mmax*mmax + 8*mmax
wa<-rep(0, nwa)
dsave<-rep(0,29)
lsave<-rep(TRUE,4)
isave<-rep(0L,44)
iwa<-rep(0L, 3*nmax)
csave<-"" # note char strings are 255 automatically


if (length(lower) == 1) lower <- rep(lower, n)
if (length(upper) == 1) upper <- rep(upper, n)

bigval <- .Machine$double.xmax/10.

for (i in 1:n) {
   if (is.finite(lower[i])) {
        if (is.finite(upper[i])) nbd[i] <- 2
        else {
           nbd[i] <- 1
           upper[i] <- bigval # to avoid call issue
             }
   } else { if (is.finite(upper[i])) {
              nbd[i] <- 3
              lower[i] <- -bigval
            } else {
              nbd[i] <- 0 
              upper[i] <- bigval
              lower[i] <- -bigval
                 }
   }
}
## cat("nbd:")
## print(nbd)


##     We start the iteration by initializing task.
## 
      itask <- 2L # START
      task <- tasklist[itask]
      f <- .Machine$double.xmax / 100
      g <- rep(f, n)

##        ------- the beginning of the loop ----------
icsave <- 0 # to make sure defined

## 111  continue ##  top of loop
repeat {
##     This is the call to the L-BFGS-B code.
      if (ctrl$trace >= 2) {
       cat("Before call, f=",f,"  task number ",itask," ")
       print(task)
      }
      result <- .Fortran('lbfgsb3', n = as.integer(n),m = as.integer(m),
                   x = as.double(prm), l = as.double(lower), u = as.double(upper),
                   nbd = as.integer(nbd), f = as.double(f), g = as.double(g),
                   factr = as.double(factr), pgtol = as.double(pgtol),
                   wa = as.double(wa), iwa = as.integer(iwa), 
                   itask = as.integer(itask),
                   iprint = as.integer(iprint),
                   icsave = as.integer(icsave), lsave=as.logical(lsave), 
                   isave=as.integer(isave), dsave=as.double(dsave))
      itask <- result$itask
      icsave <- result$icsave
      prm <- result$x
##      cat("in lbfgsb3 parameter results:")
##      print(prm)
      g <- result$g
      iwa <- result$iwa
      wa <- result$wa
      nbd <- result$nbd
      lsave <- result$lsave
      isave <- result$isave
      dsave <- result$dsave
      if (ctrl$trace > 2) {
      cat("returned from lbfgsb3\n")
      cat("returned itask is ",itask,"\n")
      task <- tasklist[itask]
      cat("changed task to ", task,"\n")
##      task<-readline("continue")
      }

      if  (itask %in% c(4L, 20L, 21L) ) {
         if (ctrl$trace >= 2) {
          cat("computing f and g at prm=")
          print(prm)
         }
##        Compute function value f for the sample problem.
         f <- fn(prm, ...)
##        Compute gradient g for the sample problem.
         if (is.null(gr)) {
             g <- grad(fn, prm, ...)
         } else {
             g <- gr(prm, ...)
         }
         if (ctrl$trace > 0) {
            cat("At iteration ", isave[34]," f =",f)
            if (ctrl$trace > 1) {
               cat("max(abs(g))=",max(abs(g)))
            }
            cat("\n")
         }
      } else {
        if (itask == 1L )  { # NEW_X
##          tmp <- readline("Continue") # eventually remove this
 		##     If task is neither FG nor NEW_X we terminate execution.
          } else break
      }
 } # end repeat
## Here build return structure
##  print(result) ## only print for debugging
  info <- list(task = task, itask = itask, lsave = lsave, 
                icsave = icsave, dsave = dsave, isave = isave)
  ans <- list(prm = prm, f = f, g = g, info = info)
##======================= The end of driver1 ============================
} # end of lbfgsb3()