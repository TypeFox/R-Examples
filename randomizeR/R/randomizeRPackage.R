#' Randomization for Clinical Trials
#' 
#' This tool enables the user to choose a randomization procedure based on 
#' sound scientific criteria. It comprises the generation of randomization 
#' sequences as well the assessment of randomization procedures based on 
#' carefully selected criteria. Furthermore, randomizeR provides a function for 
#' the comparison of randomization procedures.
#' 
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
#' 
#' @section Acknowledgement:
#' This research is embedded in the 
#' \href{http://www.ideal.rwth-aachen.de/}{IDeAl project}, which has received 
#' funding from the European Union's Seventh Framework Programme for 
#' research, technological development and demonstration under Grant Agreement 
#' no 602552. 
#' 
#' @seealso 
#' For functionality for randomization procedures, see \code{\link{randPar}} and
#' \code{\link{genSeq}}.
#' For the criteria for the assessment of randomization procedures, see
#' \code{\link{issues}}. 
#' For the assessment and comparison of randomization procedures, see 
#' \code{\link{assess}} and \code{\link{compare}}.
#' 
#' @docType package
#' @name randomizeR-package
#' @aliases randomizeR
#' @author David Schindler \email{dschindler@@ukaachen.de}, Diane Uschner \email{duschner@@ukaachen.de}, Ralf-Dieter Hilgers, Nicole Heussen
#' 
#' @import methods
#' @import ggplot2
#' @import plotrix
#' @importFrom stats dpois pt qpois qt rbinom rnorm t.test
#' @importFrom utils capture.output packageVersion sessionInfo write.table
#' @importFrom graphics abline axis box lines plot.new plot.window title
NULL
