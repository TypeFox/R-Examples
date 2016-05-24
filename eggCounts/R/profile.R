# Package relevant stuff...


minrgamma <- .Machine$double.eps  #

.verboselevel <- 0

".onAttach" <- function (lib, pkg) {
    unlockBinding(".verboselevel", asNamespace("eggCounts"))
}
