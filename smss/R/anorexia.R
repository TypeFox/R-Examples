##' @name anorexia
##'
##' @title Anorexia Study data
##' @description Weights of Anorexic girls, before and after receiving
##' one of three possible therapies. Thanks to Prof. Brian Everitt,
##' Institute of Psychiatry, London, for supplying these data.
##'
##' @docType data
##'
##' @format \Sexpr[stage=build,results=rd]{data(anorexia); smss:::describe_df(anorexia)}
##' \describe{
##'    \item{\code{subj}}{Subject ID}
##'    \item{\code{therapy}}{Therapy type. b = cognitive behavioural, f = family therapy, or c = control.}
##'    \item{\code{before}}{Weight before treatment}
##'    \item{\code{after}}{Weight after treatment}
##' }
##' @examples
##' data(anorexia)
##' summary(anorexia)
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
NULL
