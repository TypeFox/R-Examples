#'Power analysis is used in the estimation of sample sizes for experimental designs. Most programs and R packages will only output
#'the highest recommended sample size to the user. Often the user input can be complicated and computing multiple power analyses
#'for different treatment comparisons can be time consuming. This package simplifies the user input and allows the user to
#'view all of the sample size recommendations or just the ones they want to see. Currently, one-way ANOVAs \code{\link{n.oneway}}
#'and factorial ANOVAs \code{\link{n.multiway}} are supported. The effect size utilized by the functions is eta-squared which
#'is equivalent to percentage variance. It is used in the input for all of the functions so that the user may use one standard effect
#'size for all of their calculations. The calculations used to calculate the recommended sample sizes are from the 'pwr'
#'package. Future updates are planned to add more experimental designs.
#'
#'
#'
#'
#'\tabular{ll}{
#'Package: \tab easypower\cr
#'Type: \tab Package\cr
#'Version: \tab 1.0.1\cr
#'Date: \tab 2015-11-04\cr
#'License: \tab GPL (>=3)\cr
#'}
#'
#'
#'
#'
#'
#'
#'
#'@name easypower
#'@aliases easypower
#'@docType package
#'@title Sample Size Calculations Using Power Analysis
#'@author
#'\tabular{ll}{
#'Author: \tab Aaron McGarvey \email{amcgarve@@mail.uoguelph.ca}\cr
#'Maintainer: \tab Aaron McGarvey \email{amcgarve@@mail.uoguelph.ca}
#'}
#'
#'
#'
#'@importFrom "utils" "combn"
NULL
