aov.ancovaplot <- function(object, warn=TRUE) {
  object.call <- object$call
  object.call[[1]] <- as.name("aov")
  ocg <- object$call$groups
  if (!is.null(ocg)) {
    if (is.name(ocg)) {
      if (is.name(object$call[[2]][[3]])) {
        if (warn)
          warning('"groups" is specified in the call of the trellis object.\nThe formula has only one term on the rhs.\n"groups" is not included in the aov.\nPlease do the aov manually if you wish groups to be in the analysis.', call.=FALSE)
      }
      else
      object.call[[2]][[3]] <- as.call(list(as.name('+'), ocg, object.call[[2]][[3]]))
    }
    else
      warning('"groups" is specified in the call of the trellis object and "groups" is not a name.\n"groups" is not included in the aov.\nPlease do the aov manually if you wish groups to be in the analysis.', call.=FALSE)
  }
  object.call <- object.call[names(object.call) %in% c("", "data")]
  eval(object.call)
}

anova.ancovaplot <- function(object, ...) {
  anova(aov.ancovaplot(object, ...))
}

aovStatement <- function(object, ...)
  UseMethod("aovStatement")

aovStatement.ancovaplot <- function(object, ...) {
  tmp <- aov.ancovaplot(object, ...)$call
  names(tmp)[[2]] <- ""
  tmp
}

aovStatementAndAnova <- function(object, ...)
  UseMethod("aovStatementAndAnova")

aovStatementAndAnova.ancovaplot <- function(object, ...) {
  tmp <- aovStatement.ancovaplot(object, ...)
  cat(paste("> anova(", deparse(tmp), ")\n", sep=""))
  anova(object, ...)
}

model.tables.ancovaplot <- function(x, ...) {
  model.tables(aov.ancovaplot(x), ...)
}
