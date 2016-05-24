## ----setup, echo=FALSE, include=FALSE------------------------------------
library(PKNCA)

## ----eval=FALSE----------------------------------------------------------
#  library(PKNCA)
#  
#  ## Load the PK concentration data
#  d.conc <- read.csv("concentration.csv", stringsAsFactors=FALSE)
#  ## Load the dosing data
#  d.dose <- read.csv("dose.csv", stringsAsFactors=FALSE)
#  
#  ## Create a concentration object specifying the concentration, time,
#  ## study, and subject columns.  (Note that any number of grouping
#  ## levels is supporting; you are not restricted to this list.)
#  my.conc <- PKNCAconc(d.conc,
#                       concentration~time|study+subject)
#  ## Create a dosing object specifying the dose, time, study, and
#  ## subject columns.  (Note that the grouping factors should be a
#  ## subset of the grouping factors for concentration, and the columns
#  ## must have the same names between concentration and dose objects.)
#  my.dose <- PKNCAdose(d.dose,
#                       dose~time|study+subject)
#  ## Combine the concentration and dosing information both to
#  ## automatically define the intervals for NCA calculation and provide
#  ## doses for calculations requiring dose.
#  my.data <- PKNCAdata(my.conc, my.dose)
#  
#  ## Calculate the NCA parameters
#  my.results <- pk.nca(my.data)
#  
#  ## Summarize the results
#  summary(my.results)

## ------------------------------------------------------------------------
PKNCA.options()

## ----eval=FALSE----------------------------------------------------------
#  PKNCA.options(default=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  PKNCA.set.summary(name="cmax",
#                    point=business.geomean,
#                    spread=business.geocv,
#                    rounding=list(signif=3))

## ----eval=FALSE----------------------------------------------------------
#  PKNCA.set.summary(name="tmax",
#                    point=business.median,
#                    spread=business.range,
#                    rounding=list(round=2))

## ----eval=FALSE----------------------------------------------------------
#  ## Create a concentration object specifying the concentration, time,
#  ## study, and subject columns.  (Note that any number of grouping
#  ## levels is supporting; you are not restricted to this list.)
#  my.conc <- PKNCAconc(d.conc,
#                       concentration~time|study+subject/analyte)
#  ## Create a dosing object specifying the dose, time, study, and
#  ## subject columns.  (Note that the grouping factors should be a
#  ## subset of the grouping factors for concentration, and the columns
#  ## must have the same names between concentration and dose objects.)
#  my.dose <- PKNCAdose(d.dose,
#                       dose~time|study+subject)

## ------------------------------------------------------------------------
intervals <-
  data.frame(start=0, end=c(24, Inf),
             cmax=c(FALSE, TRUE),
             tmax=c(FALSE, TRUE),
             auclast=TRUE,
             aucinf=c(FALSE, TRUE))

## ----asis=TRUE, echo=FALSE-----------------------------------------------
knitr::kable(PKNCA.options()$single.dose.aucs)

## ------------------------------------------------------------------------
## find.tau can work when all doses have the same interval...
dose.times <- seq(0, 168, by=24)
print(dose.times)
PKNCA::find.tau(dose.times)

## or when the doses have mixed intervals (10 and 24 hours).
dose.times <- sort(c(seq(0, 168, by=24),
                     seq(10, 178, by=24)))
print(dose.times)
PKNCA::find.tau(dose.times)

## ----eval=FALSE----------------------------------------------------------
#  my.intervals <-
#    data.frame(start=0, end=c(24, Inf),
#               cmax=c(FALSE, TRUE),
#               tmax=c(FALSE, TRUE),
#               auclast=TRUE,
#               aucinf=c(FALSE, TRUE))
#  my.data <- PKNCAdata(my.conc, my.dose,
#                       intervals=my.intervals)

## ----eval=FALSE----------------------------------------------------------
#  my.data <- PKNCAdata(my.conc, my.dose)
#  my.intervals <- my.data$intervals
#  my.intervals$aucinf[1] <- TRUE
#  mydata$intervals <- my.intervals

