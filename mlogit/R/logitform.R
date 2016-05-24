has.intercept <- function(object, ...) {
  UseMethod("has.intercept")
}

has.intercept.default <- function(object, ...) {
  has.intercept(formula(object), ...)
}

has.intercept.formula <- function(object, ...) {
  attr(terms(object), "intercept") == 1L
}

has.intercept.Formula <- function(object, rhs = NULL, ...) {
  ## NOTE: return a logical vector of the necessary length
  ## (which might be > 1)
  if(is.null(rhs)) rhs <- 1:length(attr(object, "rhs"))
  sapply(rhs, function(x) has.intercept(formula(object, lhs = 0, rhs = x)))
}

## pFormula:
## methods : formula, model.frame, model.matrix, pmodel.response

mFormula <- function(object){
  UseMethod("mFormula")
}

is.mFormula <- function(object)
  inherits(object, "mFormula")

mFormula.formula <- function(object){
  if (!inherits(object, "Formula")) object <- Formula(object)
  class(object) <- c("mFormula", "Formula", "formula")
  object
}

mFormula <- function(object){
  stopifnot(inherits(object, "formula"))
  if (!inherits(object, "Formula")) object <- Formula(object)
  if (!inherits(object, "mFormula"))
    class(object) <- c("mFormula", class(object))
  object
}

as.Formula.mFormula <- function(x, ...){
  class(x) <- class(x)[-1]
  x
}

model.frame.mFormula <- function(formula, data, ..., lhs = NULL, rhs = NULL){
  if (is.null(rhs)) rhs <- 1:(length(formula)[2])
  if (is.null(lhs)) lhs <- ifelse(length(formula)[1]>0, 1, 0)
  index <- attr(data, "index")
##   weights <- match.call()$weights
##   print(weights)
##   if (!is.null(weights)){
##     print(class(formula))
##     formula <- as.Formula(formula(formula), paste("~", weights))
##     print(formula)
##   }
##   else{
##     formula <- as.Formula(formula)
##   }
##   mf <- model.frame(formula, as.data.frame(data))
  mf <- model.frame(as.Formula(formula), as.data.frame(data), ..., rhs = rhs)
  index <- index[rownames(mf),]
  index <- data.frame(lapply(index, function(x) x[drop = TRUE]), row.names = rownames(index))
  structure(mf,
            choice = attr(data, "choice"),
            index = index,
            class = c("mlogit.data", class(mf)))
}

has.intercept.mFormula <- function(object, ...){
  attr(object, "class") <- "Formula"
  has.int <- has.intercept(object)
  ifelse(length(has.int) > 1, has.int[2], has.int[1])
}

model.matrix.mFormula <- function(object, data, ...){
  K <- length(data)
  omitlines <- attr(na.omit(data), "na.action")
  index <- attr(data, "index")
  alt <- index[["alt"]]
  chid <- index[["chid"]]
  data$alt <- alt
  resp.name <- as.character(attr(object, "lhs")[[1]])
  # keep track of the existence of an intercept
  has.int <- has.intercept(object)
  if (has.int) intercept.char <- "alt" else intercept.char <- NULL
  
  ## for ind.spec : remove any 0 or 1 or -1 in the formula and get the
  ## list of the variables
  if (length(object)[2] > 1){
    ind.spec <- formula(object, rhs = 2, lhs = 0)
    if (!has.int) ind.spec <- update(ind.spec, ~ . + 1)
    ind.spec <- update(ind.spec, ~ .)
    ind.spec.char <- as.character(ind.spec)[2]
    if (ind.spec.char == "1") ind.spec.char <- ind.spec.var <- NULL
    else{
#      ind.spec.var <- attr(terms(ind.spec), "term.labels") the
      # following lines extract the effects and not the variable
      # names, useful for factors
      ind.spec.var <- colnames(model.matrix(update(ind.spec, ~.+1), data))[-1]
      ind.spec.char <- paste("(", ind.spec.char, "):alt", sep="")
    }
  }
  else ind.spec <- ind.spec.char <- ind.spec.var <- NULL

  # alternative specific variables
  alt.spec <- formula(object, rhs = 1, lhs = 0)
  alt.spec <- update(update(alt.spec, ~ . + 1), ~ .)
  alt.spec.char <- as.character(alt.spec)[2]
  if (alt.spec.char == "1") als.spec <- alt.spec.char <- NULL

  # specific coefficient for alternative specific variables
  if (length(object)[2] == 3){
    coef.spec <- formula(object, rhs = 3, lhs = 0)
    coef.spec <- update(update(coef.spec, ~ . + 1), ~ .)
    coef.spec.char <- as.character(coef.spec)[2]
    if (!is.null(coef.spec.char)) coef.spec.char <- paste("(", coef.spec.char, "):alt", sep="")
  }
  else coef.spec <- coef.spec.char <- NULL
  form.char <- paste(c(intercept.char, alt.spec.char,
                       ind.spec.char, coef.spec.char),
                     collapse = "+")
  formula <- as.formula(paste(resp.name, " ~ ", form.char))
  X <- model.matrix(formula, data)[, -1, drop = F]
  lev1 <- levels(alt)[1]
  lev1 <- paste("alt", lev1, sep = "")
  toremove <- unlist(lapply(as.list(ind.spec.var), function(x) paste(lev1, x, sep = ":")))
  revtoremove <- unlist(lapply(as.list(ind.spec.var), function(x) paste(x, lev1, sep = ":")))
  toremove <- colnames(X) %in% c(toremove, revtoremove)
  X <- X[, !toremove, drop = FALSE]
  # I comment the following line which seems as best to be useless
  # X[omitlines, ] <- NA

  # the following lines suppress the mentions to 'alt' in the names of
  # the effects and add a mention to '(intercept)'
  namesX <- colnames(X)
  for (i in 1:length(namesX)) namesX[i] <- sub('alt', '', namesX[i])
  z <- match(levels(alt), namesX)
  namesX[na.omit(z)] <- paste(levels(alt)[!is.na(z)], '(intercept)', sep=":")
  colnames(X) <- namesX
  X
}
