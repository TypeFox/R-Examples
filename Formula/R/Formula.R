Formula <- function(object) {

  stopifnot(inherits(object, "formula"))

  object_split <- split_formula(object)

  structure(object, lhs = object_split$lhs, rhs = object_split$rhs,
    class = c("Formula", "formula"))
}

as.Formula <- function(x, ...) UseMethod("as.Formula")

as.Formula.default <- function(x, ..., env = parent.frame()) Formula(as.formula(x, env = env))

as.Formula.Formula <- function(x, ...) x

as.Formula.formula <- function(x, ..., env) {

  ## preserve original environment
  if(missing(env)) env <- environment(x)

  ## combine all arguments to formula list
  x <- c(list(x), list(...))
  x <- lapply(x, as.formula)
  
  ## split all 
  x_split <- lapply(x, split_formula)
  x_lhs <- do.call("c", lapply(x_split, "[[", "lhs"))
  x_rhs <- do.call("c", lapply(x_split, "[[", "rhs"))

  ## recombine
  x_all <- paste_formula(x_lhs, x_rhs)
  
  ## create formula
  ## (we have everything to do this by hand, but for encapsulating code
  ## call Formula() again...which splits again)
  rval <- Formula(x_all)

  ## re-attach original environment
  environment(rval) <- env
  return(rval)
}

is.Formula <- function(object) inherits(object, "Formula")

formula.Formula <- function(x, lhs = NULL, rhs = NULL, collapse = FALSE,
  update = FALSE, drop = TRUE, ...)
{
  ## available parts
  lpart <- 1L:length(attr(x, "lhs"))
  rpart <- 1L:length(attr(x, "rhs"))

  ## default: keep all parts
  lhs <- if(is.null(lhs)) lpart else lpart[lhs]
  rhs <- if(is.null(rhs)) rpart else rpart[rhs]
  if(any(is.na(lhs))) {
    lhs <- as.vector(na.omit(lhs))
    if(length(lhs) < 1L) lhs <- 0L
    warning("subscript out of bounds, not all 'lhs' available")
  }
  if(any(is.na(rhs))) {
    rhs <- as.vector(na.omit(rhs))
    if(length(rhs) < 1L) rhs <- 0L
    warning("subscript out of bounds, not all 'rhs' available")
  }  

  ## collapse: keep parts separated by "|" or collapse with "+"
  collapse <- rep(as.logical(collapse), length.out = 2)

  rval <- paste_formula(attr(x, "lhs")[lhs], attr(x, "rhs")[rhs],
    lsep = ifelse(collapse[1L], "+", "|"),
    rsep = ifelse(collapse[2L], "+", "|"))

  ## omit potentially redundant terms
  if(all(collapse) & update) rval <- update(rval, if(length(rval) > 2) . ~ . else ~ .)

  ## reconvert to Formula if desired
  if(!drop) rval <- Formula(rval)

  ## re-attach original environment
  environment(rval) <- environment(x)

  return(rval)
}

terms.Formula <- function(x, ..., lhs = NULL, rhs = NULL, dot = "separate")
{
  ## simplify to standard formula
  form <- simplify_to_formula(x, lhs = lhs, rhs = rhs)

  ## if necessary try to expand/update/simplify formula parts with dot
  if(has_dot(form)) {
    x_orig <- x
    dot <- match.arg(dot, c("separate", "sequential"))

    ## lhs and rhs calls
    ll <- formula(x, rhs = 0L, collapse = TRUE)[[2L]]
    rr <- attr(x, "rhs")

    ## update and simplify again
    for(i in seq_along(rr)) {
      if(dot == "sequential" && i > 1L) ll <- c_formula(ll, rr[[i - 1L]], sep = "+")
      fi <- paste_formula(ll, rr[[i]]) #probably better than:# paste_formula(NULL, c_formula(rr[[i]], ll, sep = "-"))
      rr[[i]] <- update(formula(terms(fi, ...)), . ~ .)[[3L]]
    }
    attr(x, "rhs") <- rr
    form <- simplify_to_formula(x, lhs = lhs, rhs = rhs)

    ## call traditional terms()
    mt <- terms(form, ...)

    ## store updating for future reference (e.g., in model.part)
    attr(mt, "Formula_with_dot") <- x_orig
    attr(mt, "Formula_without_dot") <- x
    attr(mt, "dot") <- dot
  } else {
    ## call traditional terms()
    mt <- terms(form, ...)
  }
  
  return(mt)
}

model.frame.Formula <- function(formula, data = NULL, ..., lhs = NULL, rhs = NULL, dot = "separate")
{
  model.frame(terms(formula, lhs = lhs, rhs = rhs, data = data, dot = dot), data = data, ...)
}

model.matrix.Formula <- function(object, data = environment(object), ..., lhs = NULL, rhs = 1, dot = "separate")
{
  form <- formula(object, lhs = lhs, rhs = rhs, collapse = c(FALSE, TRUE))
  mt <- delete.response(terms(form, data = data, dot = dot))
  model.matrix(mt, data = data, ...)
}

## as model.response() is not generic, we do this:
model.part <- function(object, ...)
  UseMethod("model.part")

model.part.formula <- function(formula, data, ..., drop = FALSE) {
  formula <- Formula(formula)
  NextMethod()
}

model.part.Formula <- function(object, data, lhs = 0, rhs = 0, drop = FALSE, terms = FALSE, dot = NULL, ...) {

  ## *hs = NULL: keep all parts
  if(is.null(lhs)) lhs <- 1L:length(attr(object, "lhs"))
  if(is.null(rhs)) rhs <- 1L:length(attr(object, "rhs"))

  if(isTRUE(all.equal(as.numeric(lhs), rep(0, length(lhs)))) &
     isTRUE(all.equal(as.numeric(rhs), rep(0, length(rhs)))))
    stop("Either some 'lhs' or 'rhs' has to be selected.")

  if(is.null(dot)) {
    if(is.null(attr(attr(data, "terms"), "dot"))) {
      dot <- "separate"
    } else {
      dot <- attr(attr(data, "terms"), "dot")
    }
  } else {
    dot <- match.arg(dot, c("separate", "sequential"))
  }

  ##
  if(has_dot(object) &&
     !is.null(attr(data, "terms")) &&
     all(c("Formula_with_dot", "Formula_without_dot", "dot") %in% names(attributes(attr(data, "terms")))) &&
     dot == attr(attr(data, "terms"), "dot") &&
     simplify_to_formula(object, lhs = lhs, rhs = rhs) == simplify_to_formula(attr(attr(data, "terms"), "Formula_with_dot"), lhs = lhs, rhs = rhs)
  ) {
    object <- attr(attr(data, "terms"), "Formula_without_dot")
  }

  ## construct auxiliary terms object
  mt <- terms(object, lhs = lhs, rhs = rhs, dot = dot, data = data)

  ## subset model frame
  ix <- attr(mt, "variables")[-1L]
  if(is.null(ix)) {
    ix <- 0
  } else {
    ix <- sapply(ix, deparse)
    if(!all(ix %in% names(data))) stop(
      paste("'data' does not seem to be an appropriate 'model.frame':",
      paste(paste("'", ix[!(ix %in% names(data))], "'", sep = ""), collapse = ", "),
      "not found")
    )
  }
  rval <- data[, ix, drop = drop]
  if(!is.data.frame(rval)) names(rval) <- rownames(data)
  if(terms) attr(rval, "terms") <- mt
  return(rval)
}

update.Formula <- function(object, new,...) {

  new <- Formula(new)
  
  ## extract all building blocks
  o_lhs <- attr(object, "lhs")
  o_rhs <- attr(object, "rhs")
  n_lhs <- attr(new, "lhs")
  n_rhs <- attr(new, "rhs")
  lhs <- rep(list(NULL), length = max(length(o_lhs), length(n_lhs)))
  rhs <- rep(list(NULL), length = max(length(o_rhs), length(n_rhs)))

  ## convenience function for updating components
  update_components <- function(x, y) {
    xf <- yf <- ~ .
    xf[[2L]] <- x
    yf[[2L]] <- y
    update(xf, yf)[[2L]]
  }
    
  if(length(lhs) > 0L) for(i in 1L:length(lhs)) {
    lhs[[i]] <- if(length(o_lhs) < i) n_lhs[[i]]
      else if(length(n_lhs) < i) o_lhs[[i]]
      else update_components(o_lhs[[i]], n_lhs[[i]])
  }

  if(length(rhs) > 0L) for(i in 1L:length(rhs)) {
    rhs[[i]] <- if(length(o_rhs) < i) n_rhs[[i]]
      else if(length(n_rhs) < i) o_rhs[[i]]
      else update_components(o_rhs[[i]], n_rhs[[i]])
  }

  ## recombine
  rval <- paste_formula(lhs, rhs)
  
  ## create formula
  ## (we have everything to do this by hand, but for encapsulating code
  ## call Formula() again...which splits again)
  rval <- Formula(rval)  
  
  ## preserve original environment
  environment(rval) <- environment(object)
  
  return(rval)
}

length.Formula <- function(x) {
  ## NOTE: return length of both sides, not only rhs
  c(length(attr(x, "lhs")), length(attr(x, "rhs")))
}

print.Formula <- function(x, ...) {
  ## we could avoid calling formula() by computing on the internal
  ## structure attr(x, "rhs") <- attr(x, "lhs") <- NULL
  ## but this is probably cleaner...
  print(formula(x))
  invisible(x)
}

all.equal.Formula <- function(target, current, ...) {
  rval <- NULL
  
  if(length(target)[1L] != length(current)[1L]) {
    rval <- c(rval, paste("Length mismatch: target, current differ in number of LHS parts: ",
      length(target)[1L], ", ", length(current)[1L], sep = ""))
  } else if(!isTRUE(all.equal(attr(target, "lhs"), attr(current, "lhs")))) {
    rval <- c(rval, "Formula mismatch: LHS formulas differ in contents")
  }

  if(length(target)[2L] != length(current)[2L]) {
    rval <- c(rval, paste("Length mismatch: target, current differ in number of RHS parts: ",
      length(target)[2L], ", ", length(current)[2L], sep = ""))
  } else if(!isTRUE(all.equal(attr(target, "rhs"), attr(current, "rhs")))) {
    rval <- c(rval, "Formula mismatch: RHS formulas differ in contents")
  }
  
  if(is.null(rval)) TRUE else rval
}

str.Formula <- function(object, ...) {
  le <- length(object)
  ls <- if(sum(le) > 2L | any(le > 1L)) "s" else ""
  writeLines(c(
    sprintf("'Formula' with %s left-hand and %s right-hand side%s: %s",
      le[1L], le[2L], ls, format(object)),
    sprintf("  ..- attr(*, \".Environment\")=%s", format(attr(object, ".Environment")))))
  invisible()
}

## convenience tools #################################################

## split formulas
split_formula <- function(f) {

  stopifnot(inherits(f, "formula"))

  rhs <- if(length(f) > 2) f[[3L]] else f[[2L]]
  lhs <- if(length(f) > 2) f[[2L]] else NULL

  extract_parts <- function(x, sep = "|") {
    if(is.null(x)) return(NULL)
    
    rval <- list()
    if(length(x) > 1L && x[[1L]] == sep) {
      while(length(x) > 1L && x[[1L]] == sep) {
        rval <- c(x[[3L]], rval)
        x <- x[[2L]]
      }
    }
    return(c(x, rval))
  }

  list(lhs = extract_parts(lhs), rhs = extract_parts(rhs))
}

## combine (parts of) formulas
c_formula <- function(f1, f2, sep = "~") {

  stopifnot(length(sep) == 1L, nchar(sep) == 1L,
    sep %in% c("~", "+", "-", "|", "&"))

  if(sep == "~") {
    rval <- . ~ .
    rval[[3L]] <- f2	
    rval[[2L]] <- f1
  } else {
    rval <- as.formula(paste(". ~ .", sep, "."))
    rval[[3L]][[3L]] <- f2
    rval[[3L]][[2L]] <- f1
    rval <- rval[[3L]]
  }

  return(rval)
}

## reassemble formulas
paste_formula <- function(lhs, rhs, lsep = "|", rsep = "|") {

  stopifnot(all(nchar(lsep) == 1L), all(lsep %in% c("+", "|", "&")))
  stopifnot(all(nchar(rsep) == 1L), all(rsep %in% c("+", "|", "&")))
  
  if(length(lhs) > 1L) lsep <- rep(lsep, length.out = length(lhs) - 1L)
  if(length(rhs) > 1L) rsep <- rep(rsep, length.out = length(rhs) - 1L)

  if(is.null(lhs)) lhs <- list()
  if(is.null(rhs)) rhs <- list()
  
  if(!is.list(lhs)) lhs <- list(lhs)
  if(!is.list(rhs)) rhs <- list(rhs)

  lval <- if(length(lhs) > 0L) lhs[[1L]] else NULL
  if(length(lhs) > 1L) {
    for(i in 2L:length(lhs)) lval <- c_formula(lval, lhs[[i]], sep = lsep[[i - 1L]])
  }
  rval <- if(length(rhs) > 0L) rhs[[1L]] else 0 ## FIXME: Is there something better?
  if(length(rhs) > 1L) {
    for(i in 2L:length(rhs)) rval <- c_formula(rval, rhs[[i]], sep = rsep[[i - 1L]])
  }

  c_formula(lval, rval, sep = "~")
}

## simplify a Formula to a formula that can be processed with
## terms/model.frame etc.
simplify_to_formula <- function(Formula, lhs = NULL, rhs = NULL) {

  ## get desired subset as formula and Formula
  form <- formula(Formula, lhs = lhs, rhs = rhs)
  Form <- Formula(form)

  ## convenience functions for checking extended features
  is_lhs_extended <- function(Formula) {
    ## check for multiple parts
    if(length(attr(Formula, "lhs")) > 1L) {
      return(TRUE)
    } else {
    ## and multiple responses
      if(length(attr(Formula, "lhs")) < 1L) return(FALSE)
      return(length(attr(terms(paste_formula(NULL,
        attr(Formula, "lhs"), rsep = "+")), "term.labels")) > 1L)
    }
  }

  is_rhs_extended <- function(Formula) {
    ## check for muliple parts
    length(attr(Formula, "rhs")) > 1L
  }

  ## simplify (if necessary)
  ext_lhs <- is_lhs_extended(Form)
  if(ext_lhs | is_rhs_extended(Form)) {
    form <- if(ext_lhs) {
      if(length(attr(Form, "rhs")) == 1L & identical(attr(Form, "rhs")[[1L]], 0)) {
	paste_formula(NULL, attr(Form, "lhs"), rsep = "+")    
      } else {
        paste_formula(NULL, c(attr(Form, "lhs"), attr(Form, "rhs")), rsep = "+")
      }
    } else {
      paste_formula(attr(Form, "lhs"), attr(Form, "rhs"), rsep = "+")	 
    }
  }

  ## re-attach original environment and return
  environment(form) <- environment(Formula)
  return(form)
}

## check whether formula has a dot (FIXME: can other problems than just '.' occur?)
has_dot <- function(formula) inherits(try(terms(formula), silent = TRUE), "try-error")
