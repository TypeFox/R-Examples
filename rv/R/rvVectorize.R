

rvVectorize <- function (FUN, vectorize.args = arg.names, SIMPLIFY = FALSE, USE.NAMES = TRUE, SAMPLESIZE=NULL) {
  arg.names <- as.list(formals(FUN))
  arg.names[["..."]] <- NULL
  arg.names <- names(arg.names)
  vectorize.args <- as.character(vectorize.args)
  if (!length(vectorize.args)) {
      return(FUN)
  }
  if (!all(vectorize.args %in% arg.names)) {
      stop("must specify formal argument names to vectorize")
  }
  .FUNV <- function() {
    args <- lapply(as.list(match.call())[-1], eval, envir=parent.frame())
    names <- if (is.null(names(args))) {
      character(length(args))
    } else {
      names(args)
    }
    dovec <- names %in% vectorize.args
    Args <- .Primitive("c")(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]), 
      SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES, SAMPLESIZE=SAMPLESIZE)
    do.call(rvmapply, Args)
  }
  formals(.FUNV) <- formals(FUN)
  return(.FUNV)
}
