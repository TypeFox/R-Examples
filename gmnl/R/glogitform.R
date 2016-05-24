## gFormula:
## methods: formula, model.frame, model.matrix, pmodel.response


#' Model Formula for Multinomial Logit Models
#' 
#' Four kind of variables are used in multinomial choice models with individual heterogeneity: alternative
#' specific and individual specific variables; variables for the mean of the random parameters (deterministic taste variations), and variables
#' for the scale function. \code{gFormula} deals with this type of models using suitable methods to extract
#' the elements of the model.
#' 
#' @param object a formula for the \code{gFormula} function, for the model.matrix method, a \code{gFormula} object,
#' @param formula a \code{gFormula} object,
#' @param data a \code{data.frame},
#' @param lhs  see \code{\link[Formula]{Formula}},
#' @param rhs  see \code{\link[Formula]{Formula}},
#' @param Q  number of classes for the latent class model,
#' @param ... further arguments.
#' 
#' 
gFormula <- function(object){
  UseMethod("gFormula")
}

#' @rdname gFormula
#' @export
is.gFormula <- function(object){
  inherits(object, "gFormula")
}
  
gFormula.formula <- function(object){
  if (!inherits(object, "Formula")) object <- Formula(object)
  class(object) <- c("gFormula", "Formula", "formula")
  object
}

gFormula <- function(object){
  stopifnot(inherits(object, "formula"))
  if (!inherits(object, "Formula")) object <- Formula(object)
  if (!inherits(object, "gFormula"))
    class(object) <- c("gFormula", class(object))
  object
}

as.Formula.gFormula <- function(x, ...){
  class(x) <- class(x)[-1]
  x
}

#' @rdname gFormula
#' @import stats
#' @export
model.frame.gFormula <- function(formula, data, ..., lhs = NULL, rhs = NULL){
  if (is.null(rhs)) rhs <- 1:(length(formula)[2])
  if (is.null(lhs)) lhs <- ifelse(length(formula)[1] > 0, 1, 0)
  mf <- model.frame(as.Formula(formula), as.data.frame(data), ..., rhs = rhs)
  index <- attr(data, "index")
  index <- index[rownames(mf), ]
  index <- data.frame(lapply(index , function(x) x[drop = TRUE]), rownames(index))
  structure(mf,
            choice = attr(data, "choice"),
            index = index,
            class = c("mlogit.data", class(mf)))
}

#' @rdname gFormula
#' @import stats
#' @export
model.matrix.gFormula <- function(object, data, rhs = NULL, Q = NULL, ...){
  K <- length(data) # Number of parameters
  index <- attr(data, "index")
  alt <- index[["alt"]]
  chid <- index[["chid"]]
  data$alt <- alt
  resp.name <- as.character(attr(object, "lhs"))
  
  if (is.null(rhs)) {
    has.int <- has.intercept(object)
    if (has.int) intercept.char <- "alt" else intercept.char <- NULL
    
    ## for ind.spec : remove any 0 or 1 or -1 in the formula and get the
    ## list of the variables
    
    if (length(object)[2] > 1) {
      ind.spec <- formula(object, rhs = 2, lhs = 0) # individual specific variables
      if (!has.int) ind.spec <- update(ind.spec, ~. + 1)
      ind.spec <- update(ind.spec, ~.)
      ind.spec.char <- as.character(ind.spec)[2]
      if (ind.spec.char == "1") ind.spec.char <- ind.spec.var <- NULL
      else {
        ind.spec.var <- colnames(model.matrix(update(ind.spec, ~. + 1), data))[-1]
        ind.spec.char <- paste("(", ind.spec.char, "):alt", sep = "")
      }
    }
    else ind.spec <- ind.spec.char <- ind.spec.var <- NULL
    
    # alternative specific variables
    alt.spec <- formula(object, rhs = 1, lhs = 0)
    alt.spec <- update(update(alt.spec, ~. + 1), ~.)
    alt.spec.char <- as.character(alt.spec)[2]
    if (alt.spec.char == "1") als.spec <- alt.spec.char <- NULL
    
    # specific coefficient for alternative specific variables
    if (length(object)[2] > 2) {
      coef.spec <- formula(object, rhs = 3, lhs = 0)
      coef.spec <- update(update(coef.spec, ~. + 1), ~.)
      coef.spec.char <- as.character(coef.spec)[2]
      if (!is.null(coef.spec.char)) coef.spec.char <- paste("(", coef.spec.char, "):alt", sep = "")
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
    
    # the following lines suppress the mentions to 'alt' in the names of
    # the effects and add a mention to '(intercept)'
    namesX <- colnames(X)
    for (i in 1:length(namesX)) namesX[i] <- sub('alt', '', namesX[i])
    z <- match(levels(alt), namesX)
    namesX[na.omit(z)] <- paste(levels(alt)[!is.na(z)], '(intercept)', sep = ":")
    colnames(X) <- namesX
  } else {
    if (rhs < 4) stop("rhs should be greater than 3")
    id <- index[["id"]]
    if (is.null(id)) {
      indata <- data[!duplicated(chid), ]
    } else {
      indata <- data[!duplicated(id), ]
    }
    if (is.null(Q)) {
      ind.var <- formula(object, rhs = rhs, lhs = 0)
      X <- model.matrix(ind.var, indata)
    } else {
      cldata   <- indata[rep(seq_len(nrow(indata)), each = Q), ] # expand data Q times
      if (is.null(id)) {
        chid.un <- unique(chid)
        class   <- factor(rep(1:Q, length(chid.un)))
      } else {
        id.un <- unique(id)
        class   <- factor(rep(1:Q, length(id.un)))
      }
      cldata   <- cbind(cldata, class, row.names = NULL)
      cldata   <- mlogit::mlogit.data(cldata, choice = resp.name, alt.var = "class",  shape = "long")
      index <- attr(cldata, "index")
      alt <- index[["alt"]]
      class.var <- formula(object, rhs = rhs, lhs = 0) 
      has.int <- has.intercept(class.var)
      if (has.int) intercept.char <- "factor(class)" else intercept.char <- NULL 
      if (!has.int) {
        class.var <- update(class.var, ~. + 1)
        class.var<- update(class.var, ~.)
        class.var.char <- as.character(class.var)[2]
        if (class.var.char == "1") class.var.char <- class.var.var <- NULL
      } else {
         has.xclass <- as.character(class.var)[2]
        if (has.xclass == "1") {
          class.var.char <- NULL
        } else {
          class.var.var <- colnames(model.matrix(update(class.var, ~. + 1), cldata))[-1]
          class.var.char <- paste("(", class.var.var, "):class", sep = "")
        } 
      }
      form.char <- paste(c(intercept.char, class.var.char), collapse = "+")
      formula <- as.formula(paste(resp.name, " ~ ", form.char))
      X <- model.matrix(formula, cldata)[, -1, drop = F]
      
      lev1 <- levels(class)[1]
      lev1 <- paste("class", lev1, sep = "")
      if (has.xclass != "1") {
        toremove <- unlist(lapply(as.list(class.var.var), function(x) paste(lev1, x, sep = ":")))
        revtoremove <- unlist(lapply(as.list(class.var.var), function(x) paste(x, lev1, sep = ":")))
        toremove <- colnames(X) %in% c(toremove, revtoremove)
        X <- X[, !toremove, drop = FALSE]
      } 
      namesX <- colnames(X)
      for (i in 1:length(namesX)) namesX[i] <- sub('factor', '', namesX[i])
      colnames(X) <- namesX
      attr(X, "alt") <- alt
    }
  }
  X
}

## has.intercept
has.intercept <- function(object, ...) {
  UseMethod("has.intercept")
}

#' @import stats
has.intercept.default <- function(object, ...) {
  has.intercept(formula(object), ...)
}

#' @import stats
has.intercept.formula <- function(object, ...) {
  attr(terms(object), "intercept") == 1L
}

#' @import stats
has.intercept.Formula <- function(object, rhs = NULL, ...) {
  if (is.null(rhs)) rhs <- 1:length(attr(object, "rhs"))
  sapply(rhs, function(x) has.intercept(formula(object, lhs = 0, rhs = x)))
}

has.intercept.gFormula <- function(object, ...){
  attr(object, "class") <- "Formula"
  has.int <- has.intercept(object)
  ifelse(length(has.int) > 1, has.int[2], has.int[1])
}

## has.othervar
has.othervar <- function(object, ...) {
  UseMethod("has.othervar")
}

has.othervar.default <- function(object, ...) {
  has.othervar(object, ...)
}

has.othervar.Formula <- function(object, ...) {
  therhs <- attr(object, "rhs")
  therhs
}

has.othervar.gFormula <- function(object, rhs, ...){
  attr(object, "class") <- "Formula"
  therhs <- has.othervar(object)
  len <- length(object)[2]
  if (len < rhs || therhs[[rhs]] == "0") FALSE else TRUE
}

