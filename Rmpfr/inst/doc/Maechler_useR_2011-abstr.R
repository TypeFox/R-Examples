### R code from vignette source 'Maechler_useR_2011-abstr.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        width = 75)
Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")
stopifnot(require("Rmpfr"))


###################################################
### code chunk number 2: ex-exp
###################################################
options(digits = 17)# to print to full "standard R" precision
.N <- function(.) mpfr(., precBits = 200)

exp(   1 )
exp(.N(1))


###################################################
### code chunk number 3: nice-but-does-not-fit-on-1-page (eval = FALSE)
###################################################
## choose    ( 200, 99:100 )
## chooseMpfr( 200, 99:100 )


