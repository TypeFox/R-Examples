#' @name cdtamodel-class
#' @title Class cdtamodel
#' @description A cdtamodel class in the CopulaDTA package.
#' @docType class
#' @slot copula Copula function, 'fgm', 'gauss', 'c90', '270', or 'frank'.
#' @slot modelcode Character with the model code as returned by the model function
#' @slot modelargs List containing control parameters for the prior distributions
#' @family cdta
#' @seealso \link{cdtamodel}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}

setClass(Class="cdtamodel",
          representation=representation(
              copula = 'character',
              modelcode = 'character',
              modelargs = "list"))

#' @name cdtafit-class
#' @title Class cdtafit
#' @description A cdtafit class in the CopulaDTA package.
#' @docType class
#' @slot data A data-frame with no missing values containg TP, TN, FP, FN, 'SID' and co-varaiables(if necessary).
#' @slot SID A string indicating the name of the column with the study identifier.
#' @slot copula Copula function, 'fgm', 'gauss', 'c90', '270', or 'frank'.
#' @slot modelargs List containing control parameters for the prior distributions.
#' @slot fit An object of class stanfit returned by the function sampling.
#' @family cdta
#' @seealso \link{fit}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
#' @importClassesFrom rstan stanmodel stanfit


setClass(Class="cdtafit",
         representation=representation(
            data='data.frame',
            SID = 'character',
            copula = 'character',
            modelargs = "list",
            fit='stanfit'))
