## model.frame.mvr: Extract or generate the model frame from a `mvr' object.
## It is simply a slightly modified `model.frame.lm'.
model.frame.groc <- function(formula, ...)
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    fcall[[1]] <- as.name("mvr")
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) env <- parent.frame()
    eval(fcall, env, parent.frame())
  }
  else formula$model
}

## model.matrix.mvr: Extract the model matrix from an `mvr' object.
## It is a modified version of model.matrix.lm.
model.matrix.groc <- function(object, ...)
{
  if (n_match <- match("x", names(object), 0))
    object[[n_match]]
  else {
    data <- model.frame(object, ...)
    mm <- NextMethod("model.matrix", data = data)
    mm <- delete.intercept(mm) # Deletes any intercept coloumn
    ## model.matrix.default prepends the term name to the colnames of
    ## matrices.  If there is only one predictor term, and the
    ## corresponding matrix has colnames, remove the prepended term name:
    mt <- terms(object)
    if (length(attr(mt, "term.labels")) == 1 &&
        !is.null(colnames(data[[attr(mt, "term.labels")]])))
      colnames(mm) <- sub(attr(mt, "term.labels"), "", colnames(mm))
    return(mm)
  }
}


## delete.intercept: utilitiy function that deletes the response coloumn from
## a model matrix, and adjusts the "assign" attribute:
delete.intercept <- function(mm) {
  ## Save the attributes prior to removing the intercept coloumn:
  saveattr <- attributes(mm)
  ## Find the intercept coloumn:
  intercept <- which(saveattr$assign == 0)
  ## Return if there was no intercept coloumn:
  if (!length(intercept)) return(mm)
  ## Remove the intercept coloumn:
  mm <- mm[,-intercept, drop=FALSE]
  ## Update the attributes with the new dimensions:
  saveattr$dim <- dim(mm)
  saveattr$dimnames <- dimnames(mm)
  ## Remove the assignment of the intercept from the attributes:
  saveattr$assign <- saveattr$assign[-intercept]
  ## Restore the (modified) attributes:
  attributes(mm) <- saveattr
  ## Return the model matrix:
  mm
}

