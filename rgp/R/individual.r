## individual.R
##   - Functions for analysing GP individuals
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Patrick Koch, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Functions for analysing GP individuals
##'
##' \code{inputVariablesOfIndividual} returns a list of input variables in \code{inset} that are
##' used by the GP individual \code{ind}.
##'
##' @param ind A GP individual, represented as a R function.
##' @param inset A set of input variables.
##'
##' @rdname individualAnalysis
##' @export
inputVariablesOfIndividual <- function(ind, inset)
  as.list(intersect(FlattenExpression(body(ind)), inset))
