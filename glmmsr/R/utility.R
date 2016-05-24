has_re <- function(formula) {
  length(lme4::findbars(formula[[length(formula)]])) > 0L
}

find_dim_in_dir <- function(x, i) {
  d <- length(dim(x))
  if(d == 0 && i == 1) {
    # have a vector
    return(length(x))
  } else if(d < i) {
    stop(paste0("Array does not have dimension ", i), call. = FALSE)
  } else {
    return(dim(x)[i])
  }
}

subset_dim <- function(x, i, indices){
  d <- length(dim(x))
  if(d == 0 && i == 1) {
    # have a vector
    return(x[indices])
  }
  if(d < i) {
    stop(paste0("Array does not have dimension ", i), call. = FALSE)
  }
  if(i > 1L) {
    commas_before <- paste(rep(",", times = i - 1), collapse = " ")
  } else {
    commas_before <- character(0)
  }
  if(i < d) {
    commas_after <- paste(rep(",", times = d - i), collapse = " ")
  } else {
    commas_after <- character(0)
  }
  expr_text <- paste("x[", commas_before, "indices", commas_after,
                     ", drop = FALSE]")
  eval(parse(text = expr_text))
}


# copied from terms.formula
# stop dropping of brackets around (1 | group) terms
fixFormulaObject <- function(object) {
  Terms <- terms(object)
  tmp <- attr(Terms, "term.labels")
  ind <- grep("|", tmp, fixed = TRUE)
  if (length(ind))
    tmp[ind] <- paste("(", tmp[ind], ")")
  if (length(ind <- attr(Terms, "offset"))) {
    tmp2 <- as.character(attr(Terms, "variables"))[-1L]
    tmp <- c(tmp, tmp2[ind])
  }
  rhs <- if (length(tmp))
    paste(tmp, collapse = " + ")
  else "1"
  if (!attr(terms(object), "intercept"))
    rhs <- paste(rhs, "- 1")
  if (length(form <- formula(object)) > 2L) {
    res <- formula(paste("lhs ~", rhs))
    res[[2L]] <- form[[2L]]
    res
  }
  else formula(paste("~", rhs))
}


find_pairs <- function(set) {
  if(length(set) > 1L) {
    return(utils::combn(set, m = 2, simplify = FALSE))
  } else {
    return(NULL)
  }
}

check_weights <- function(weights) {
  tol <- .Machine$double.eps^0.5
  if(length(weights) > 0) {
    weights_are_integers <- all(abs(weights - round(weights)) < tol)
    if(!weights_are_integers)
      stop("Cannot currently handle non-integer weights", call. = FALSE)
  }
}

check_modfr_SR <- function(modfr) {
  family <- modfr$family$family
  link <- modfr$family$link
  if(family != "binomial")
    stop("Only binomial family currently implemented for sequential reduction approximation",
         call. = FALSE)
  if(!link %in% c("logit", "probit"))
    stop("Only logit and probit links currently implemented for sequential reduction approximation",
         call. = FALSE)
}
