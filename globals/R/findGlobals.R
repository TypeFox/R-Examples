## This function is equivalent to:
##    fun <- asFunction(expr, envir=envir, ...)
##    codetools::findGlobals(fun, merge=TRUE)
## but we expand it here to make it more explicit
## what is going on.
#' @importFrom codetools makeUsageCollector findLocalsList walkCode
findGlobals_conservative <- function(expr, envir, ...) {
  objs <- character()

  enter <- function(type, v, e, w) {
    objs <<- c(objs, v)
  }

  ## From codetools::findGlobals():
  fun <- asFunction(expr, envir=envir, ...)
  # codetools::collectUsage(fun, enterGlobal=enter)

  ## The latter becomes equivalent to (after cleanup):
  w <- makeUsageCollector(fun, enterGlobal=enter, name="<anonymous>")
  w$env <- new.env(hash=TRUE, parent=w$env)
  locals <- findLocalsList(list(expr))
  for (name in locals) assign(name, value=TRUE, envir=w$env)
  walkCode(expr, w)

  unique(objs)
}


#' @importFrom codetools makeUsageCollector walkCode
findGlobals_liberal <- function(expr, envir, ...) {
  objs <- character()

  enter <- function(type, v, e, w) {
    objs <<- c(objs, v)
  }

  fun <- asFunction(expr, envir=envir, ...)

  w <- makeUsageCollector(fun, enterGlobal=enter, name="<anonymous>")
  walkCode(expr, w)

  unique(objs)
}


#' @importFrom codetools makeUsageCollector walkCode
findGlobals_ordered <- function(expr, envir, ...) {
  class <- name <- character()

  enterLocal <- function(type, v, e, w) {
    class <<- c(class, "local")
    name <<- c(name, v)
  }

  enterGlobal <- function(type, v, e, w) {
    class <<- c(class, "global")
    name <<- c(name, v)
  }

  fun <- asFunction(expr, envir=envir, ...)

  w <- makeUsageCollector(fun, name="<anonymous>",
                          enterLocal=enterLocal, enterGlobal=enterGlobal)
  walkCode(expr, w)

  ## Drop duplicated names
  dups <- duplicated(name)
  class <- class[!dups]
  name <- name[!dups]

  unique(name[class == "global"])
}


#' @export
findGlobals <- function(expr, envir=parent.frame(), ..., tweak=NULL, dotdotdot=c("warning", "error", "return", "ignore"), method=c("ordered", "conservative", "liberal"), substitute=FALSE, unlist=TRUE) {
  method <- match.arg(method)
  dotdotdot <- match.arg(dotdotdot)

  if (substitute) expr <- substitute(expr)

  if (is.list(expr)) {
    globals <- lapply(expr, FUN=findGlobals, envir=envir, ..., tweak=tweak, dotdotdot=dotdotdot, substitute=FALSE, unlist=FALSE)
    if (unlist) {
      needsDotdotdot <- FALSE
      for (kk in seq_along(globals)) {
        s <- globals[[kk]]
        n <- length(s)
        if (identical(s[n], "...")) {
          needsDotdotdot <- TRUE
          s <- s[-n]
          globals[[kk]] <- s
        }
      }
      globals <- unlist(globals, use.names=TRUE)
      globals <- sort(unique(globals))
      if (needsDotdotdot) globals <- c(globals, "...")
    }
    return(globals)
  }

  if (is.function(tweak)) expr <- tweak(expr)

  if (method == "ordered") {
    findGlobalsT <- findGlobals_ordered
  } else if (method == "conservative") {
    findGlobalsT <- findGlobals_conservative
  } else if (method == "liberal") {
    findGlobalsT <- findGlobals_liberal
  }

  ## Is there a need for global '...' variables?
  needsDotdotdot <- FALSE
  globals <- withCallingHandlers({
    oopts <- options(warn=0L)
    on.exit(options(oopts))
    findGlobalsT(expr, envir=envir)
  }, warning=function(w) {
    ## Warned about '...'?
    pattern <- "... may be used in an incorrect context"
    if (grepl(pattern, w$message, fixed=TRUE)) {
      needsDotdotdot <<- TRUE
      if (dotdotdot == "return") {
        ## Consume / muffle warning
        invokeRestart("muffleWarning")
      } else if (dotdotdot == "ignore") {
        needsDotdotdot <<- FALSE
        ## Consume / muffle warning
        invokeRestart("muffleWarning")
      } else if (dotdotdot == "error") {
        e <- simpleError(w$message, w$call)
        stop(e)
      }
    }
  })

  if (needsDotdotdot) globals <- c(globals, "...")

  globals
}
