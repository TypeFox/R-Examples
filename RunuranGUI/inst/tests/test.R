#############################################################################
##
##  Test Runuran GUI
##
#############################################################################

## -- Constants -------------------------------------------------------------

## if TRUE, then run tests in development mode
DEVEL <- FALSE
DEVEL <- TRUE

## if TRUE, then all tests are performed
ALL <- TRUE
##ALL <- FALSE

## if TRUE, then dialog boxes are non-modal (use galert instead of gmessage)
## (works in developement mode only!)
TESTING <- TRUE
##TESTING <- FALSE

## time unit (speed) for running tests
time.unit <- 1

## -- Load required libraries -----------------------------------------------

if (DEVEL) {
  cat("******************************\n",
      "*  Run in development mode.  *\n",
      "******************************\n\n", sep="")
  
  ## load libraries
  require(Runuran)
  require(rvgtest)
  require(gWidgets)
  options("guiToolkit"="RGtk2")
  ##options("guiToolkit"="tcltk")

  ## for development we simple 'source' the the code
  source("../../R/unuran_gui.R", local=FALSE)
  source("../../R/distributions.R", local=FALSE)
  source("../../R/methods.R", local=FALSE)
  source("../../R/stage1.R", local=FALSE)
  source("../../R/stage2.R", local=FALSE)
  source("../../R/stage3.R", local=FALSE)
  source("../../R/hist.R", local=FALSE)
  source("../../R/aerror.R", local=FALSE)
  source("../../R/common.R", local=FALSE)
  source("../../R/emulate.R", local=FALSE)

  RunuranGUI.TESTING <- TESTING

} else {
  cat("****************************\n",
      "*  Run in installed mode.  *\n",
      "****************************\n\n", sep="")

  ## load installed library
  require(RunuranGUI)
  ## we need non-exported functions
  emul.stage1 <- RunuranGUI:::emul.stage1
  emul.stage2 <- RunuranGUI:::emul.stage2

  if (isTRUE(TESTING))
    warning("Testing mode does not work")
}

## -- Global variable -------------------------------------------------------

## time delay between each triggering events
delay <- time.unit

## start timer
time.start <- proc.time()

#############################################################################


## -- Stage 1: Start, restart and cancel ------------------------------------

if (ALL) {
  w <- unuran.gui()
  emul.stage1(w, delay=delay, distr=list(type="discr", howdef="user-defined"))
  emul.stage1(w, delay=delay, button="new")
  emul.stage1(w, delay=delay, button="cancel")
}

## -- Stage 1: Error: missing data ------------------------------------------

if (ALL) {
  w <- unuran.gui()
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont", howdef="built-in"))
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="discr", howdef="built-in"))
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont", howdef="built-in", cont="Normal"),
              method=list(type="Select method"))
  emul.stage1(w, delay=delay, button="cancel")
}

## -- Stage 2: Start, restart and cancel ------------------------------------

if (ALL) {
  w <- unuran.gui()
  emul.stage1(w, delay=delay, distr=list(type="cont",cont="Normal"))
  emul.stage1(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="new")
  emul.stage1(w, delay=delay, distr=list(type="cont",cont="Normal"))
  emul.stage1(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="cancel")
}

## -- Stage 2: Parameters of built-in distributions -------------------------

if (ALL) {
  w <- unuran.gui()
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Automatic"))
  ## fail
  emul.stage2(w, delay=delay, button="ok",
              distr=list(sd=-1))
  emul.stage2(w, delay=delay, button="ok",
              distr=list(sd=1, VAR3=1, VAR4=0)) ## lower and upper boundary

  emul.stage2(w, delay=delay, button="cancel")
}

## -- Stage 2: Generation method for continuous built-in distributions ------

if (ALL) {
  w <- unuran.gui()

  ## -- Automatic ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Automatic"))
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="new")

  ## -- Inversion ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Inversion"))
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="new")

  ## -- Rejection ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Rejection"))
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="new")

  ## -- ARS ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Select method", cont="ARS"))
  ## default:
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="new")

  ## fail
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Gamma"),
              method=list(type="Select method", cont="ARS"))
  emul.stage2(w, delay=delay, button="ok",
              distr=list(shape=0.5))
  emul.stage2(w, delay=delay, button="new")

  ## -- ITDR ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Gamma"),
              method=list(type="Select method", cont="ITDR"))
  ## default:
  emul.stage2(w, delay=delay, button="ok",
              distr=list(shape=0.5))
  emul.stage2(w, delay=delay, button="new")

  ## fail
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Select method", cont="ITDR"))
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="new")

  ## -- PINV ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Select method", cont="PINV"))
  ## default:
  emul.stage2(w, delay=delay, button="ok")
  ## invalid:
  emul.stage2(w, delay=delay, button="ok",
              method=list(uresolution="junk", smooth=TRUE))
  ## ok:
  emul.stage2(w, delay=delay, button="ok",
              method=list(uresolution=1e-12, smooth=FALSE))
  emul.stage2(w, delay=delay, button="new")

  ## -- SROU ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Select method", cont="SROU"))
  ## default:
  emul.stage2(w, delay=delay, button="ok")
  ## invalid:
  emul.stage2(w, delay=delay, button="ok",
              method=list(r="junk"))
  ## ok:
  emul.stage2(w, delay=delay, button="ok",
              method=list(r=2))
  emul.stage2(w, delay=delay, button="new")

  ## fail
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Gamma"),
              method=list(type="Select method", cont="SROU"))
  emul.stage2(w, delay=delay, button="ok",
              distr=list(shape=0.5))
  emul.stage2(w, delay=delay, button="new")

  ## -- TDR ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Select method", cont="TDR"))
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="new")

  ## fail
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont",cont="Gamma"),
              method=list(type="Select method", cont="TDR"))
  emul.stage2(w, delay=delay, button="ok",
              distr=list(shape=0.5))
  emul.stage2(w, delay=delay, button="new")


  ## -- stop ---
  emul.stage1(w, delay=delay, button="cancel")
}

## -- Stage 2: Generation method for continuous user-defined distributions --

if (ALL) {
  w <- unuran.gui()

  ## -- ARS ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont", howdef="user-defined"),
              method=list(type="Select method", cont="ARS"))
  ## fail
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="-x"))        ## log-density
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR2="-1"))        ## derivative
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR4=0))           ## lower bound
  emul.stage2(w, delay=delay, button="new")

  ## -- ITDR ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont", howdef="user-defined"),
              method=list(type="Select method", cont="ITDR"))

  ## fail
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="-0.5*log(x)-x"))   ## log-density
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR2="-0.5/x-1"))        ## derivative
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR3=TRUE))              ## is log-density
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR4=0))                 ## pole
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR5=0))                 ## lower bound
  emul.stage2(w, delay=delay, button="new")

  ## -- PINV ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont", howdef="user-defined"),
              method=list(type="Select method", cont="PINV"))

  ## fail
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="exp(-x)", VAR3=FALSE))   ## density
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR2="exp(-x)"))               ## CDF
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR4=0))                       ## mode
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR5=0, VAR6=2))               ## domain
  emul.stage2(w, delay=delay, button="new")

  ## -- SROU ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont", howdef="user-defined"),
              method=list(type="Select method", cont="SROU"))

  ## fail
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="-x", VAR2=TRUE))   ## log-density
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR3=0, VAR4=1, VAR5=0, VAR6=Inf)) ## mode, area, domain
  emul.stage2(w, delay=delay, button="new")

  ## -- TDR ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="cont", howdef="user-defined"),
              method=list(type="Select method", cont="TDR"))

  ## fail
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="-x", VAR2="-1", VAR3=TRUE))   ## log-density
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR4=0, VAR5=Inf))  ## domain
  emul.stage2(w, delay=delay, button="new")

  ## -- stop ---
  emul.stage1(w, delay=delay, button="cancel")

}

## -- Stage 2: Generation method for built-in discrete distributions --------

if (ALL) {
  w <- unuran.gui()

  ## -- DAU ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="discr", discr="Binomial"),
              method=list(type="Select method", discr="DAU"))
  ## fail
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="ok",
              distr=list(size=100, prob=0.3))
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR4=100))                 ## domain
  emul.stage2(w, delay=delay, button="new")

  ## -- DGT ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="discr", discr="Binomial"),
              method=list(type="Select method", discr="DGT"))
  ## fail
  emul.stage2(w, delay=delay, button="ok")
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(size=100, prob=0.3, VAR4=100))
  emul.stage2(w, delay=delay, button="new")

  ## -- DARI ---
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="discr", discr="Binomial"),
              method=list(type="Select method", discr="DARI"))
  ## fail
  emul.stage2(w, delay=delay, button="ok")
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(size=100, prob=0.3, VAR4=100))
  emul.stage2(w, delay=delay, button="new")

  ## -- stop ---
  emul.stage1(w, delay=delay, button="cancel")

}

## -- Stage 2: Generation method for built-in discrete distributions --------

if (ALL) {

  ## -- DARI ---
  w <- unuran.gui()
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="discr", howdef="user-defined"),
              method=list(type="Select method", discr="DARI"))
  ## fail
  emul.stage2(w, delay=delay, button="ok")
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="dbinom(x,100,0.3)",lb=0,ub=100))   ## pmf
  emul.stage2(w, delay=delay, button="ok",
              distr=list(mode="30"))
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(sum="1"))
  emul.stage2(w, delay=delay, button="cancel")
  
  ## -- DAU ---
  w <- unuran.gui()
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="discr", howdef="user-defined"),
              method=list(type="Select method", discr="DAU"))
  ## fail
  emul.stage2(w, delay=delay, button="ok")
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="c(1,2,3,4,5)"))  ## pv
  emul.stage2(w, delay=delay, button="cancel")

  ## -- DGT ---
  w <- unuran.gui()
  emul.stage1(w, delay=delay, button="ok",
              distr=list(type="discr", howdef="user-defined"),
              method=list(type="Select method", discr="DGT"))
  ## fail
  emul.stage2(w, delay=delay, button="ok")
  ## ok
  emul.stage2(w, delay=delay, button="ok",
              distr=list(VAR1="c(1,2,3,4,5)"))  ## pv
  emul.stage2(w, delay=delay, button="cancel")

}

## -- Stage 2: Perform ------------------------------------------------------

if (ALL) {
  w <- unuran.gui()
  emul.stage1(w, delay=delay, distr=list(cont="Normal"), button="ok")

  emul.stage2(w, delay=delay, button="ok",
              action=list(
                distr=c(TRUE,"mydist1"),
                gen=c(TRUE,"mygen1"),
                sample=c(TRUE,"mysample1",1234)))
  emul.stage2(w, delay=delay, button="ok",
              action=list(
                distr=c(TRUE,"mydist2"),
                gen=c(FALSE,"mygen2"),
                sample=c(FALSE,"mysample2",2345)))
  emul.stage2(w, delay=delay, button="ok",
              action=list(
                distr=c(FALSE,"mydist3"),
                gen=c(TRUE,"mygen3"),
                sample=c(FALSE,"mysample3",3456)))
  emul.stage2(w, delay=delay, button="ok",
              action=list(
                distr=c(FALSE,"mydist4"),
                gen=c(FALSE,"mygen4"),
                sample=c(TRUE,"mysample4",4567)))

  emul.stage2(w, delay=delay, button="ok",
              action=list(
                distr=c(TRUE,"my()dist2"),
                gen=c(FALSE,"mygen2"),
                sample=c(FALSE,"mysample2",2345)))
  emul.stage2(w, delay=delay, button="ok",
              action=list(
                distr=c(FALSE,"mydist3"),
                gen=c(TRUE,"my()gen3"),
                sample=c(FALSE,"mysample3",3456)))
  emul.stage2(w, delay=delay, button="ok",
              action=list(
                distr=c(FALSE,"mydist4"),
                gen=c(FALSE,"mygen4"),
                sample=c(TRUE,"my()sample4",4567)))

  emul.stage2(w, delay=delay, button="cancel")
}

#############################################################################

if (0) {
  w <- unuran.gui()

  ## -- Automatic ---
  emul.stage1(w, delay=0, button="ok",
              distr=list(type="cont",cont="Normal"),
              method=list(type="Automatic"))
  emul.stage2(w, delay=100*delay, button="cancel")
}


run.time <- (proc.time() - time.start)[3]  ## "elapsed" time
run.time

#############################################################################

#emul.stage2(w, delay=delay, action=list(notebook=TRUE))
## --------------------------------------------------------------------------

