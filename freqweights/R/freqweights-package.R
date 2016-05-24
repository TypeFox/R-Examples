

##' Working with frequency tables
##' 
##' The frequency of a particular data value is the number of times it occurs.
##' A frequency table is a table of values with their corresponding
##' frequencies. Frequency weights are integer numbers that indicate how many
##' cases each case represents. This package provides some functions to work
##' with such type of collected data. 
##' 
##' \tabular{ll}{
##' Package: \tab freqweights\cr
##' Type: \tab Package\cr
##' Version: \tab 0.1.0\cr
##' Date: \tab 2014-05-20\cr
##' License: \tab GPL 3.0\cr }
##' 
##' @name freqweights-package
##' @aliases freqweights-package freqweights
##' @docType package
##' @author Emilio Torres-Manzanera
##' 
##' Maintainer: Emilio Torres-Manzanera <torres@@uniovi.es>
##' 
##' @keywords package
##' @examples
##' tablefreq(iris)
##' lmfreq(Sepal.Length ~ Petal.Length, tablefreq(iris))
##' hclustvfreq(tablefreq(iris[,1:4]))
NULL
