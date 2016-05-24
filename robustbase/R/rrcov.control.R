##  rrcov : Scalable Robust Estimators with High Breakdown Point
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, a copy is available at
##  http://www.r-project.org/Licenses/


## "FIXME": If you change this, you must "sync" with ../man/rrcov.control.Rd
##          1) covMcd()'s         default in ./covMcd.R
##          2) ltsReg.default()'s default in ./ltsReg.R
##          3) covComed()s        default in ./comedian.R
rrcov.control <-
    function(alpha = 1/2, method = c("covMcd", "covComed", "ltsReg"),
             nsamp = 500, nmini = 300, kmini = 5,
             seed = NULL, tolSolve = 1e-14,
             scalefn = "hrv2012", maxcsteps = 200,
	     trace = FALSE, wgtFUN = "01.original", beta,
             use.correction = identical(wgtFUN, "01.original"),
             adjust = FALSE)
{
    method <- match.arg(method)
    if(missing(beta) || !is.numeric(beta))
        beta <- c("covMcd" = 0.975, "ltsReg" = 0.9875, "covComed" = 0.95)[[method]]
    list(alpha=alpha, nsamp=nsamp, nmini=as.integer(nmini), kmini=as.integer(kmini),
         seed = as.integer(seed),
	 tolSolve=tolSolve, scalefn=scalefn, maxcsteps=as.integer(maxcsteps),
         trace=trace, wgtFUN=wgtFUN, beta=beta,
	 use.correction=use.correction, adjust=adjust)
}
## allow direct fast access:
.scalefn.default <- eval(formals(rrcov.control)$scalefn)

## Only for back compatibility, as some new args did not exist pre 2013-04,
## and callers of ltsReg() / covMcd() may use a "too small"  'control' list:
getDefCtrl <- function(nm, defCtrl = rrcov.control()) {
    callerEnv <- parent.frame()
    if(is.null(get(nm, envir = callerEnv)))
	assign(nm, defCtrl[[nm]], envir=callerEnv)
}
