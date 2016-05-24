##---------------------------------------------------------------------------
## A Physiologically Based Pharmacokinetic (PBPK) model
## before trying this code, the C or FORTRAN program has to be compiled
## this can be done in R:
## system("R CMD SHLIB satres.f")
## or:
## system("R CMD SHLIB satresC.c")
## do make sure that this file is in the working directory...
## (if not, use setwd() )
##---------------------------------------------------------------------------

## We want to be able to run three kinds of dosing regimens with the same
## model:
##  - single gavage
##  - repeated gavage
##  - dietary 

library(deSolve)

wh <- menu(c("C version", "FORTRAN version"),
           graphics = TRUE,
           title = "Which language version?")
if (wh == 0) stop("User cancelled", cal. = FALSE)
DLLname <- switch(wh, "satresC", "satres")
FullDLLname <- paste(DLLname, .Platform$dynlib.ext, sep = "")
if (!file.exists(FullDLLname)) stop(paste("You need to create", FullDLLname,
                                          "using 'R CMD SHLIB", DLLname, "'"),
                                    call. = FALSE)
dyn.load(FullDLLname)

if(length(grep("intakes", search())) == 0) attach("intakes.RData")

## Dose is the Dose in mg/kg
## Doseint is NA for single dose, interval
## between doses in hours for repeated dosing, -1 to use the intake data
## MaxTime is the largest requested output time, and is calculated
## internally.

## Other parms as in satres.c
defParms <-
  c(Vc = 0.0027, Vt = 0.0545, kd = 0.00059/0.0027, ka = 0.537,
    Tm = 860.9, KT = 0.0015, kfil = 0.6830/0.0027,
    Vfil = 0.01, free = 0.02, BW = 0.025, Dose = NA,
    Doseint = NA, Qd = NA, Qfil = NA, MaxTime = NA, TDose = NA)

## initparms is called as, for example
## P <- initparms(list(Dose = 60, Doseint = 24, Vc = 0.0030))
## Gives a parameter list that the model can use, for 60 mg/kg
## every 24 hours dosing, setting Vc to 0.003 L
initparms <- function(newParms = NULL) {
  Parms <- defParms
  if (!is.null(newParms)) {
    ldots <- as.list(newParms)
    if (!all(names(ldots) %in% names(defParms))) 
      stop("illegal parameter name")
    Parms[names(ldots)] <- unlist(ldots)    
  }
  lParms <- as.list(Parms)
  Parms["Qd"] <- with(lParms, kd * Vc)
  Parms["Qfil"] <- with(lParms, kfil*Vc)
  Parms["TDose"] <- Parms["Dose"] * Parms["BW"]
  Parms
}

## newParms is a list with parameter names
initforcs <- function(Parms) {
  if (is.na(Parms["Doseint"]))
    RepDose <- matrix(c(0, Parms["MaxTime"], 0, 0), ncol = 2)
    
  else if (Parms["Doseint"] > 0) {
    Parms["TDose"] <- Parms["TDose"]/(5/3600)
    dosetimes <- seq(0, Parms["MaxTime"] - 5/3600, by = Parms["Doseint"])
    dosesoff <- dosetimes + 5/3600
    RepDose <- cbind(sort(c(dosetimes, dosesoff)),
                     rep(c(Parms["TDose"], 0), length(dosetimes)))

  } else if (Parms["Doseint"] < 0) {
    maxdays <- ceiling(Parms["MaxTime"]/24)
    dosetimes <- as.vector(outer(intakes[, "hours"], 24*(0:maxdays), "+"))
    doserates <- rep(intakes[, "Rfood.femaleB6C3F1"], (maxdays + 1)) * Parms["TDose"]
    RepDose <- cbind(dosetimes, doserates)
  }
  RepDose
}

## initState returns the initialized state vector.
initstate <- function(Parms){
  if (is.na(Parms["Doseint"])) 
    structure(c(rep(0, 3), Parms["TDose"], 0, 0),
              names = c("Ccentral", "Csecond", "Cfiltrate", "Agut", "Elim",
                "AUC"))
  else 
    structure(rep(0, 6),
              names = c("Ccentral", "Csecond", "Cfiltrate", "Agut", "Elim",
                "AUC"))
}

## pfoasat runs the model.  On input,
##   - Times is a vector of time values
##     at which model results are desired.
##   - newParms is a list like the input
##     to initparms, above.
##   - method is a string giving the solution method to use
##     see the documentation for deSolve::ode for details
##     there.  the elipsis (...) is for additional arguments
##     to the odesolver (see ode and the individual methods
##     for details).
## The return value is a matrix of values.  Column 1 is the
## time vector, Columns 2 - 5 are the concentrations in
## compartments 1 - 4 (just before dosing, in the case of repeated
## dosing).
##
## Example: to match the 7 and 17 day 20 mg/kg repeated dosing
## using lsode:
## out <- pfoasat(24 * c(0, 7, 17), newParms = list(Dose = 20, Doseint = 24))
## when finished, you can unload the dll with
## dyn.unload("satres")
pfoasat <- function(Times, newParms, method = "lsode", ...){
  
  if ("MaxTime" %in% names(newParms)) 
    newParms["MaxTime"] <- max(Times)
  else 
    newParms <- c(newParms, MaxTime = max(Times))

  Parmsout <- initparms(newParms)
  Forcings <- initforcs(Parmsout)
  y <- initstate(Parmsout)

  ode(y, Times, "derivs", parms = Parmsout, method = method,
      dllname = DLLname, initfunc = "initmod", forcings = Forcings,
      initforc = "initforc", fcontrol = list(method = "constant"),
      nout = 1, outnames = "Total", ...)
}

## -------------------------------------------------------------------
## Simulate a range of doses, both be repeated gavage and an equivalent
## dose via the diet.  Plot the time course for 1 and 500 mg/kg/day,
## and the total dose-response.
Doses <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
nperhour <- 6 ## for smooth plotting
ndays <- 30 ## follow for ndays

outs <- vector("list", length = 2*length(Doses))
dim(outs) <- c(length(Doses), 2)
rownames(outs) <- as.character(Doses)

for (i in seq(along = Doses)) {
  
outs[[i, 1]] <- list(Dose = Doses[i],
                    out = as.data.frame(pfoasat(seq(0, 24*ndays, by = 1/nperhour),
                      newParms = list(Dose = Doses[i], Doseint = 24),
                      hmax = 0.001))
                      )

outs[[i, 2]] <- list(Dose = Doses[i],
                    out = as.data.frame(pfoasat(seq(0, 24*ndays, by = 1/nperhour),
                      newParms = list(Dose = Doses[i], Doseint = -1),
                      hmax = 0.4))
                    )
}

## Plot 1 and 500 mg/kg/day doses, to see the contrast
par(mfrow = c(1, 2), las = 1, bty = "l", mar = c(5, 4, 0, 1))

## ------------------------ Central compartment
ylim = c(0, 500)
plot(Ccentral ~ I(time/24), data = outs[["1", 1]]$out, type = "l", ylim = ylim,
     xlab = "Days in Study",
     ylab = "Conc. PFOA in Central Cmpt.",
     sub = "A: 1 mg/kg/day")
lines(Ccentral ~ I(time/24), data = outs[["1", 2]]$out, lty = "44")
legend("right",legend = c("Daily gavage", "Feed"), lty = c("solid", "44"),
       bty = "n")
plot(Ccentral ~ I(time/24), data = outs[["500", 1]]$out, type = "l", ylim = ylim,
     xlab = "Days in Study",
     ylab = "Conc. PFOA in Central Cmpt.",
     sub = "B: 500 mg/kg/day")
lines(Ccentral ~ I(time/24), data = outs[["500", 2]]$out, lty = "44")

## Force a pause after this figure
tmp <- readline(prompt = "press <Enter> to continue ... ")

## now, the curve relating external dose to internal dose-metric.

## Function to extract the dose-metric
## z is a dataframe like the ones we've made here
## We compute the average daily peak concentration in
## the central compartment and the daily average AUC in the
## central compartment.
dosemetric <- function(z) {
  ## drop the first time (0)
  z <- z[-1,]
  ## split the data on day:
  day <- ceiling(z$time/24)
  dailypeaks <- tapply(z$Ccentral, day, function(x) max(x))
  dailyaucs <- tapply(z$AUC, day, function(x) (x[length(x)] - x[1]))/24
  c(avgpeak = mean(dailypeaks), avgauc = mean(dailyaucs))
}

## Create a matrix to hold the doses
DoseMets <- matrix(nrow = length(Doses), ncol = 4,
                   dimnames = list(rownames(outs),
                     c("gavage.peak", "gavage.auc", "diet.peak", "diet.auc")))
for (dose in rownames(outs)) {
  DoseMets[dose, c("gavage.peak", "gavage.auc")] <- dosemetric(outs[[dose, 1]]$out)
  DoseMets[dose, c("diet.peak", "diet.auc")] <- dosemetric(outs[[dose, 2]]$out)
}

DoseMets <- as.data.frame(cbind(Doses, DoseMets))

par(mfrow = c(1, 1), bty = "l", las = 1, mar = c(4, 4, 0, 0))
plot(gavage.peak ~ Doses, DoseMets,
     ylim = range(DoseMets[, 2:5]),
     xlab = "Administered Dose (mg/kg/day)",
     ylab = "Dose Metric", log = "xy", pch = 1)
zz <- spline(log(DoseMets$Doses), log(DoseMets$gavage.peak))
lines(exp(zz[[1]]), exp(zz[[2]]))
points(gavage.auc ~ Doses, DoseMets, pch = 20)
zz <- spline(log(DoseMets$Doses), log(DoseMets$gavage.auc))
lines(exp(zz[[1]]), exp(zz[[2]]), lty = "33")
points(diet.peak ~ Doses, DoseMets, pch = 1, col = "blue")
zz <- spline(log(DoseMets$Doses), log(DoseMets$diet.peak))
lines(exp(zz[[1]]), exp(zz[[2]]), col = "blue")
points(diet.auc ~ Doses, DoseMets, pch = 20, col = "blue")
zz <- spline(log(DoseMets$Doses), log(DoseMets$diet.auc))
lines(exp(zz[[1]]), exp(zz[[2]]), lty = "33", col = "blue")
legend("topleft", legend = c("gavage peak", "gavage AUC", "diet peak", "diet auc"),
       pch = c(1, 20, 1, 20), lty = c("solid", "33", "solid", "33"),
       col = c("black", "black", "blue", "blue"),
       bty = "n")

## unload the DLL
dyn.unload(FullDLLname)
