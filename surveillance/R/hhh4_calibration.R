################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### calibrationTest() for "hhh4" fits
###
### Copyright (C) 2015 Sebastian Meyer
### $Revision: 1428 $
### $Date: 2015-07-16 12:14:52 +0200 (Don, 16. Jul 2015) $
################################################################################

calibrationTest.hhh4 <- function (x,
                                  subset = x$control$subset,
                                  units = seq_len(x$nUnit),
                                  ...)
{
    ## perform the calibration test in the specified subset
    res <- calibrationTest.default(
        x = x$stsObj@observed[subset, units, drop = FALSE],
        mu = x$fitted.values[match(subset, x$control$subset), units, drop = FALSE],
        size = psi2size.hhh4(x, subset, units),
        ...)

    ## change "data.name" to be the name of the supplied model
    res$data.name <- deparse(substitute(x))
    res
}

calibrationTest.oneStepAhead <- function (x, ...)
{
    ## perform the calibration test
    res <- calibrationTest.default(
        x = x$observed,
        mu = x$pred,
        size = psi2size.oneStepAhead(x),
        ...)

    ## change "data.name" to be the name of the supplied "oneStepAhead" object
    res$data.name <- deparse(substitute(x))
    res

}
