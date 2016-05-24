## rgp.R
##   - RGP package documentation for Roxygen
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' RGP is a simple yet flexible modular Genetic Programming system for the R
##' environment. The system implements classical untyped tree-based genetic
##' programming as well as more advanced variants including, for example,
##' strongly typed genetic programming and Pareto genetic programming.
##' 
##' @docType package
##' @name rgp-package
##' @title The RGP package
##' @author Oliver Flasch \email{oliver.flasch@@fh-koeln.de}, Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}, Thomas Bartz-Beielstein \email{thomas.bartz-beielstein@@fh-koeln.de}, Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}, Joerg Stork \email{joerg.stork@@fh-koeln.de}
##' @keywords package
##' @useDynLib rgp

NULL

## Package startup functions used to initialize static objects...
.onLoad <- function(libname, pkgname) {
  # initialize types for standard GP functions...
  "rnorm" %::% (list(st("numeric"), st("numeric"), st("numeric")) %->% st("numeric"))
  "runif" %::% (list(st("numeric"), st("numeric"), st("numeric")) %->% st("numeric"))
  "+" %::% (list(st("numeric"), st("numeric")) %->% st("numeric"))
  "-" %::% (list(st("numeric"), st("numeric")) %->% st("numeric"))
  "*" %::% (list(st("numeric"), st("numeric")) %->% st("numeric"))
  "/" %::% (list(st("numeric"), st("numeric")) %->% st("numeric"))
  "sqrt" %::% (list(st("numeric")) %->% st("numeric"))
  "exp" %::% (list(st("numeric")) %->% st("numeric"))
  "log" %::% (list(st("numeric")) %->% st("numeric"))
  "sin" %::% (list(st("numeric")) %->% st("numeric"))
  "cos" %::% (list(st("numeric")) %->% st("numeric"))
  "tan" %::% (list(st("numeric")) %->% st("numeric"))

  # initialize standard GP function and constant sets...
  arithmeticFunctionSet <<- functionSet("+", "-", "*", "/")
  expLogFunctionSet <<- functionSet("sqrt", "exp", "log")
  trigonometricFunctionSet <<- functionSet("sin", "cos", "tan")
  mathFunctionSet <<- c(arithmeticFunctionSet, expLogFunctionSet, trigonometricFunctionSet)

  numericConstantSet <<- constantFactorySet(function() rnorm(1, 0, 1))
}

.onAttach <- function(libname, pkgname) {
  # show startup message
  packageStartupMessage("*** RGP version ", (sessionInfo())$otherPkg$rgp$Version, " initialized successfully.\n",
                        "    Type 'help(package=\"rgp\")' to bring up the RGP help pages,\n",
                        "    or type 'vignette(\"rgp_introduction\")' to show RGP's package vignette.\n",
                        "    Type 'symbolicRegressionUi()' to bring up the symbolic regression UI if\n",
                        "    the optional package 'rgpui' is installed.")
}

## ... .
