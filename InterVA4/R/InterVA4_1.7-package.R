

#' Translation list of COD codes
#' 
#' This is the translation of COD abbreviation codes into their corresponding
#' full names.
#' 
#' 
#' @name causetext
#' @docType data
#' @format A data frame with the translation of COD codes to their names on 68
#' CODs (both the version of COD only and COD with group code).
#' @keywords datasets
#' @examples
#' 
#' data(causetext)
#' 
NULL





#' Perform InterVA4 algorithm and provide graphical summarization of COD
#' distribution.
#' 
#' Computes individual cause of death and population cause-specific mortality
#' fractions using the InterVA4 algorithm. Provides a simple graphical
#' representation of the result.
#' 
#' To get the most up-to-date version of the package, as well as the past
#' versions, please check the github repository at:
#' \url{https://github.com/richardli/InterVA4/}
#' 
#' \tabular{ll}{ Package: \tab InterVA4\cr Type: \tab Package\cr Version: \tab
#' 1.6\cr Date: \tab 2015-08-29\cr License: \tab GPL-2\cr }
#' 
#' @name InterVA4-package
#' @aliases InterVA4 InterVA4-package
#' @docType package
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @references http://www.interva.net/
#' @keywords InterVA
#' @examples
#' 
#' data(SampleInput)
#' sample.output <- InterVA(SampleInput, HIV = "h", Malaria = "v", directory = "VA test", 
#'     filename = "VA_result", output = "extended", append = FALSE)
#' 
NULL





#' Conditional probability of InterVA4
#' 
#' This is the table of conditional probabilities of symptoms given CODs. The
#' values are from InterVA-4.1.
#' 
#' 
#' @name probbase
#' @docType data
#' @format A data frame with 246 observations on 81 variables. Each observation
#' is the conditional probability.
#' @keywords datasets
#' @examples
#' 
#' data(probbase)
#' 
NULL





#' 10 records of Sample Input
#' 
#' This is a dataset consisting of 10 arbitrary sample input deaths in the
#' acceptable format of InterVA4. Any data that needs to be analyzed by this
#' package should be in the same format. The orders of the input fields must
#' not be changed.
#' 
#' 
#' @name SampleInput
#' @docType data
#' @format 10 arbitrary input records.
#' @keywords datasets
#' @examples
#' 
#' data(SampleInput)
#' 
NULL



