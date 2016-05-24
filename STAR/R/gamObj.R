gamObj <- function(object,
                   ...) {

  UseMethod("gamObj")
}

gamObj.gampsth <- function(object,...) {

  evalq(PoissonF, envir=environment(object$lambdaFct))

}

gamObj.gamlockedTrain <- function(object, ...) {
  object$gamFit
}
