##'--------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2016-01-10
##'
##' Description:
##' New generics for exported methods
##'
##' Overview:
##' cvloss:             cross-validates the loss of a model
##'                     regarding a tuning parameter. See
##'                     'tvcm-cv.R'.
##' extract:            extracts features of a fitted model.
##'                     See 'tvcm-methods.R'.
##' neglogLik2:         extracts the -2*log-likelihood loss of
##'                     a fitted model. See also 'olmm-methods.R'
##'                     and 'tvcm-methods.R'
##' neglogLik2.default: neglogLik2 function which should
##'                     work for most model classes
##' oobloss:            estimates out-of-bag loss.
##' otsplot:            ordinal time series plot.
##' prunepath:          extracts the path of pruning a fitted
##'                     tree.
##' ranefCov:           extracts covariance matrix of random
##'                     effect variance parameters.
##' splitpath:          extracts the splitting path of the
##'                     growing process of a tree.
##'
##' Last modifications:
##' 2014-09-07: added 'prunepath' generic
##' 2014-07-17: Update the descriptions
##'--------------------------------------------------------- #

cvloss <- function(object, ...) UseMethod("cvloss")

extract <- function(object, ...) UseMethod("extract")

fixef.glm <- function(object, ...) coef(object)

neglogLik2 <- function(object, ...) UseMethod("neglogLik2")

neglogLik2.default <- function(object, ...)
    return(-as.numeric(2 * logLik(object)))

oobloss <- function(object, ...) UseMethod("oobloss")

otsplot <- function(x, ...) UseMethod("otsplot")

prunepath <- function(tree, ...) UseMethod("prunepath")

ranefCov <- function(object, ...) UseMethod("ranefCov")

splitpath <- function(tree, ...) UseMethod("splitpath")
