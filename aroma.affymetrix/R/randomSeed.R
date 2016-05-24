###########################################################################/**
# @RdocFunction randomSeed
#
# @title "Sets and resets the .Random.seed in the global environment"
#
# \description{
#  @get "title".
# }
#
# \usage{
#  @usage randomSeed
# }
#
# \arguments{
#   \item{action}{A @character string specifying the action.}
#   \item{seed}{Random seed to be set; only for \code{action="set"}.
#     If \code{length(seed) == 1}, then \code{set.seed(seed)} is
#     used, otherwise \code{.Random.seed} is assigned the value.}
#   \item{kind}{(optional) A @character string specifying type of
#     random number generator to use, cf. @see "base::RNGkind".}
#   \item{backup}{If @TRUE, the previous (seed, kind) state is recorded
#     such that it can be reset later.}
# }
#
# \value{
#   Returns a \code{.Random.seed}.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
randomSeed <- local({
  oldSeed <- NULL
  oldKind <- NULL
  lecuyerSeed <- NULL

  genv <- globalenv()

  getSeed <- function() {
    if (exists(".Random.seed", envir=genv, inherits=FALSE)) {
      get(".Random.seed", envir=genv, inherits=FALSE)
    } else {
      NULL
    }
  }

  setSeed <- function(seed, kind=NULL, backup=TRUE) {
    force(seed)  ## FIX: Why is this needed?

    ## Set new RNG kind?
    newKind <- (!is.null(kind) && !identical(kind, RNGkind()[1L]))
    if (newKind) {
       if (backup) {
         oldSeed <<- getSeed()
         oldKind <<- RNGkind()[1L]
       }
       RNGkind(kind)  ## Sets .Random.seed
    }

    ## Reset or set seed?
    if (is.null(seed)) {
      if (exists(".Random.seed", envir=genv, inherits=FALSE)) {
        rm(list=".Random.seed", envir=genv, inherits=FALSE)
        lecuyerSeed <<- NULL
      }
    } else {
      if (backup && !newKind) oldSeed <<- getSeed()

      if (length(seed) == 1L) {
        set.seed(seed)
        lecuyerSeed <<- getSeed()
      } else {
        assign(".Random.seed", seed, envir=genv, inherits=FALSE)
        lecuyerSeed <<- seed
      }
    }
  }

  advanceSeed <- function() {
    ## Nothing to do?
    if (RNGkind()[1L] != "L'Ecuyer-CMRG") return()

    if (is.null(lecuyerSeed)) {
      stats::runif(1)
      lecuyerSeed <<- getSeed()
    }

    lecuyerSeed <<- nextRNGStream(lecuyerSeed)
    assign(".Random.seed", lecuyerSeed, envir=genv, inherits=FALSE)
  }


  function(action=c("set", "advance", "reset", "get"), seed=NULL, kind=NULL, backup=TRUE) {
    action <- match.arg(action)

    ## Record existing RNG kind (only once)
    if (is.null(oldKind)) oldKind <<- RNGkind()[1L]

    if (action == "set") {
      setSeed(seed=seed, kind=kind, backup=backup)
    } else if (action == "advance") {
      advanceSeed()
    } else if (action == "reset") {
      setSeed(seed=oldSeed, kind=oldKind, backup=FALSE)
    } else if (action == "get") {
      return(getSeed())
    }

    invisible(getSeed())
  }
}) # randomSeed()
