##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2016-01-10
##'
##' Description:
##' General utility functions for the 'vcrpart' package.
##' most of the functions are not exported and documented inline.
##'
##' Functions:
##' neglogLik2.glm:        compute the -2 times likelihood error of
##'                        a 'glm' object.
##' addEmptyChar:          add empty spaces to a character
##' appendDefArgs:         over-write default arguments
##' deparseCall:           convert a 'call' into a 'character'
##' formatMatrix:          format matrices for print functions
##' vcrpart_copy:          duplicate R objects
##' vcrpart_value_space:   extract the values space of data data.frame
##' fe, vc, re, ce, ge     special terms for formulas. Are exported.
##' vcrpart_contr.sum:     compute weighted sum contrasts
##' vcrpart_fitted:        extract model fitted values
##' vcrpart_formula_eta:   constructs a formula for linear predictors
##' vcrpart_formula_cond:  constructs a formula for contitioning variables
##' vcrpart_formula:       extracts a list of predictor formulas
##' vcrpart_formula_delEnv: deletes environments of lists of formulas
##'                         from vcrpart_formula.
##'
##' Last modifications:
##' 2015-09-02: - modified 'vcrpart_formula_eta' for implementation of
##'               gaussian mixed model.
##'             - set for baseline model as default category-specific
##'               effects. 
##' 2015-06-01: corrected bug of 'vcrpart_formula' in cases where a
##'             effect modifier has the variable name 'x'. Now I call
##'             it 'fTerm'.
##' 2015-01-20: change defaults for random effect specification
##' 2014-11-10: exported function 'contr.wsum'
##' 2014-09-08: partial substitution of 'rep' by 'rep.int'
##' 2014-09-07: added 'vcrpart_fitted' function to avoid replicated
##'             definitions of the same function.
##' 2014-09-04: removed 'adj' option in 'vc'
##' 2014-08-29: add the 'adj' argument to 'vc' to adjust the selection
##'             if the number of predictors per vc() term varies.
##' 2014-08-27: allow character vectors in 'vc' terms
##' 2014-07-31: - added to possibility to remove the intercept
##'               directly in the formula by a '-1'
##'             - now more than one varying intercepts are
##'               allowed to be specificed (the only case in which
##'               this was critical was when using 'sctest = TRUE',
##'               and this is not checked in 'tvcm' itself)
##' 2014-06-17: modified documentation style:
##' 2014-06-17: new function 'vcrpart_copy' that duplicates objects
##' 2014-06-17: removed 'vcrpart_slot' after converting 'olmm'
##'             class to a S3 class
##' 2014-06-03: modify 'tvcm_formula' to allow component-wise
##'             trees
##'
##' To do:
##' - 
##' -------------------------------------------------------- #

## --------------------------------------------------------- #
##' Add empty spaces to characters to have the same length.
##'
##' @param x       character vector.
##' @param nchar   integer. The size per string.
##' @param justify "left" or "right". The place where the spaces
##'                should be added.
##'
##' @return A character vector.
## --------------------------------------------------------- #

addEmptyChar <- function(x, nchar, justify = "left") {
  
  for (i in 1:length(x))
    if (justify == "left") {
      x[i] <- paste(x[i], paste(rep.int(" ", nchar - nchar(x[i])),
                                collapse = ""), sep = "")
    } else if (justify == "right") {
      x[i] <- paste(paste(rep.int(" ", nchar - nchar(x[i])),
                          collapse = ""), x[i], sep = "")
    }
  return(x)
}


## --------------------------------------------------------- #
##' Over-writes default- with user-specified arguments
##'
##' @param args    a list of user-specified arguments.
##' @param default a list of default arguments.
##'
##' @return A list with arguments.
## --------------------------------------------------------- #

appendDefArgs <- function(args, default) {
  
  if (is.null(args)) return(default)
  subs <- setdiff(names(default), names(args))
  if (length(subs) > 0)
    for (i in subs) args[[i]] <- default[[i]] 
  return(args)
}


## --------------------------------------------------------- #
##' Converts a call slot into a character. The function is
##' used for summary and print functions.
##'
##' @param x slot of a call, e.g., call$formula.
##'
##' @return A character.
## --------------------------------------------------------- #

deparseCall <- function(x) {
  
  if (is.null(x)) return(character())
  rval <- paste(deparse(x), collapse = "\n")
  if (grepl("structure(list(", rval, fixed = TRUE)) rval <- character()
  return(rval)
}


## --------------------------------------------------------- #
##' Format colnames etc. of a matrix for pretty prints.
##'
##' @param x a matrix
##'
##' @return A matrix.
## --------------------------------------------------------- #

formatMatrix <- function(x, ...) {
  
  rowNames <- rownames(x)
  colNames <- colnames(x)
  if (is.null(rowNames)) rowNames <- rep.int("", nrow(x))
  if (!is.null(colNames)) rowNames <- append("", rowNames)
  rowNames <- addEmptyChar(rowNames, max(nchar(rowNames)))
  rval <- format(x, ...)
  if (!is.null(colNames)) {
    for (i in 1L:ncol(rval))
      if (nchar(colNames[i]) < max(nchar(rval[i])))
        colNames[i] <- addEmptyChar(colNames[i], max(nchar(rval[i])), "right")
    for (i in 1L:ncol(rval)) 
      for (j in 1L:nrow(rval)) 
        if (nchar(rval[j, i]) < nchar(colNames[i]))
          rval[j, i] <- addEmptyChar(rval[j, i], nchar(colNames[i]), "right")
    rval <- rbind(colNames, rval)
  }
  rval <- cbind(rowNames, rval)
  rval <- apply(rval, 1, paste, collapse = "  ")
  rval <- paste(rval, "\n", sep = "")
  rval <- paste(rval, collapse = "")
  return(rval)
}


## --------------------------------------------------------- #
##' Duplicate R objects
##'
##' Duplicates an R object such that any modification on the
##' values of the copied object will not change the values
##' of the original object. This function may be used to
##' protect values of the source object when modifying the
##' copied object within a '.Call' call.
##' 
##' @param x the R object to be duplicated.
##' 
##' @return A list with arguments.
##'
##' @details used in 'tvcm_fit_loss'
## --------------------------------------------------------- #

vcrpart_copy <- function(x)
  return(.Call("vcrpart_duplicate", x, PACKAGE = "vcrpart"))


## --------------------------------------------------------- #
##' Extract the values space of data data.frame.
##'
##' @param data  a data.frame.
##' @param neval the maximum number of values to be evaluated
##'    for numeric vectors. 
##'
##' @return A list with values of the variables in 'data'
## --------------------------------------------------------- #

vcrpart_value_space <- function(data, neval = 50L) {
  
  FUN <- function(x, neval) {
    rval <- sort(unique(x))
    if (is.numeric(x) && length(rval) > neval) {
      rval <- as.double(quantile(rval, seq(0,1,length.out = neval)))
    }
    if (is.factor(x)) rval <- droplevels(rval)
    return(rval)
  }
  return(lapply(as.list(data), FUN, neval = neval))
}


fe <- function(formula, intercept = TRUE) {
  stopifnot(is.logical(intercept) | is.character(intercept))
  stopifnot(length(intercept) == 1L)
  if (is.logical(intercept))
    intercept <- ifelse(intercept, "ce", "none")
  mc <- match.call()
  if (missing(formula)) {
    formula <- if (intercept == "none") "-1" else "1"
  } else {
    formula <- deparse(mc$formula, 500L)
  }
  eta <- formula(paste("~", formula))
  return(list(eta = eta, cond = ~1, intercept = intercept, type = "fe"))
}


vc <- function(..., by, intercept = missing(by), 
               nuisance = character()) {
  stopifnot(is.logical(intercept) | is.character(intercept))
  stopifnot(is.character(nuisance))
  if (is.logical(intercept))
    intercept <- ifelse(intercept, "ce", "none")
  if (!intercept %in% c("none", "ce", "ge"))
    stop("'intercept' in 'fe' must be a logical or one of 'none', 'ce' or 'ge'.")
  mc <- match.call()
  eta <- if (missing(by)) formula(~1) else formula(paste("~", deparse(mc$by, 500L)))
  if ((intercept == "ce" & length(all.vars(eta)) == 0L) |
      (intercept == "none")) nuisance <- setdiff(nuisance, "(Intercept)")
  nuisance <- intersect(nuisance, c("(Intercept)", all.vars(eta)))
  subs <- if (is.null(names(mc))) 2:length(mc) else which(names(mc) == "")[-1L]
  if (length(subs) < 1L) stop("no effect modifiers in 'vc'.")
  cond <- "~"
  for (i in subs) {
    if (i != subs[1L]) cond <- paste(cond, "+")
    isChar <- inherits(try(eval.parent(mc[[i]]), silent = TRUE), "character")
    if (isChar) {
      cond <- paste(cond, paste(eval.parent(mc[[i]]), collapse = "+"))
    } else {
      cond <- paste(cond, deparse(mc[[i]], 500L))
    }
  }
  cond <- formula(cond)
  return(list(eta = eta, cond = cond, intercept = intercept,
              nuisance = nuisance, type = "vc"))
}


re <- function(formula, intercept = TRUE) {
  stopifnot(is.logical(intercept) | is.character(intercept))
  if (is.logical(intercept))
    intercept <- ifelse(intercept, "ge", "none")
  if (!intercept %in% c("none", "ce", "ge"))
    stop("'intercept' in 'fe' must be a logical or one of 'none', 'ce' or 'ge'")
  mc <- match.call()
  formula <- gsub(" ", "", deparse(mc$formula, 500L))
  formula <- strsplit(formula, "|", fixed = TRUE)[[1L]]
  if (length(formula) < 2L) stop("no grouping factor 'subject' in 're'.")
  if (length(formula) > 2L) stop("'formula' in 're' is missspecified.")
  cond <- formula(paste("~", formula[2]))
  if (length(all.vars(formula(paste("~", cond)))) > 1L)
    stop("maximum one grouping factor 'subject' is allowed in 're'")
  if (formula[1L] == "") formula[1L] <- "1"
  eta <- formula(paste("~", formula[1L]))
  return(list(eta = eta, cond = cond, intercept = intercept, type = "re"))
}


ce <- function(formula) {
  mc <- match.call()
  formula <- deparse(mc$formula, 500L)
  formula <- formula(paste("~", formula))
  return(attr(terms(formula, keep.order = TRUE), "term.labels"))
}


ge <- function(formula) {
  mc <- match.call()
  formula <- deparse(mc$formula, 500L)
  formula <- formula(paste("~", formula))
  return(attr(terms(formula, keep.order = TRUE), "term.labels"))
}

contr.wsum <- function(x, weights = rep.int(1.0, length(x))) {
  stopifnot(is.factor(x))
  stopifnot(is.numeric(weights))
  stopifnot(length(x) == length(weights))
  con <- NULL
  if (nlevels(x) > 1L) {
    con <- contr.sum(levels(x))
    tab <- tapply(weights, x, sum)
    con[nrow(con),] <- con[nrow(con),] * tab[-length(tab)] / tab[length(tab)]
    colnames(con) <- levels(x)[1:(nlevels(x) - 1)]
  }
  return(con)
}


## --------------------------------------------------------- #
##' Extract model fitted values
##' 
##' Extracts model fitted values. The function essentially
##' calls 'predict' and deletes 'newdata' if supplied
##'
##' @param object  a fitted object.
##' @param ...     further arguments passed to predict
##'
##' @return A matrix of fitted values.
##'
##' @details Substitutes 'fitted.fvcm', 'fitted.tvcm' and
##'    'fitted.olmm'
## --------------------------------------------------------- #

vcrpart_fitted <- function(object, ...) {
  call <- call(name = "predict", object = object)
  dotargs <- list(...)
  dotargs$newdata <- NULL
  for (arg in names(dotargs)) call[[arg]] <- list(...)[[arg]]
  return(eval(call))
}


## --------------------------------------------------------- #
##' Utility function for \code{vcrpart_formula}
##' 
##' Constructs a formula for the linear predictor.
##'
##' @param x       a formula.
##' @param fit     the fitting function, e.g., olmm() or glm()
##' @param env     environment for evaluating the formula.
##'
##' @return A formula.
##'
##' @details Used in \code{vcrpart_formula}.
## --------------------------------------------------------- #

vcrpart_formula_eta <- function(x, family, env) {

  terms <- terms(x$eta, specials = c("fe", "vc", "re", "ce", "ge"), keep.order = TRUE)
  if (length(unlist(attr(terms, "specials")[c("fe", "vc", "re")])) > 0)
    stop("'", x$type, "' term contains 'fe', 'vc' or 're' terms.")
  rval <- attr(terms, "term.labels")
  termFact <- rownames(attr(terms, "factors"))
  
  if (x$type %in% c("fe", "re", "vc") & inherits(family, "family.olmm")) {

    ## extract the terms for 'olmm' objects which may include
    ## terms width 'ge()' and 'ce'

    ## extract terms
    if (length(rval) > 0L) {

      checkOperators <- function(x) {
        x <- gsub(" ", "", x)
        return(grepl(":", x, fixed = TRUE) |
               grepl("*", x, fixed = TRUE) |
               grepl("%in%", x))
      }
      
      subsCe <- rval %in% termFact[attr(terms, "specials")$ce]
      ceTerms <-
          unlist(lapply(rval[subsCe], function(fTerm) eval(parse(text = fTerm))))
      
      if (x$type == "vc" && any(sapply(ceTerms, checkOperators)))
        stop("the ':', '*' and '%in%' operators are not allowed for the 'by' ",
             "argument in 'vc' terms.")
      
      subsGe <- rval %in% termFact[attr(terms, "specials")$ge]
      geTerms <-
          unlist(lapply(rval[subsGe], function(fTerm) eval(parse(text = fTerm)))) 
          
      if (x$type == "vc" && any(sapply(geTerms, checkOperators)))
        stop("the ':', '*' and '%in%' operators are not allowed for the 'by' ",
             "argument in 'vc' terms.")

      ## This code distinguishes between the models and their random
      ## effect specification. For the cumulative and the adjacent model,
      ## use global effects, and for the cumulative model use category
      ## specific effects
      if (family$family %in% c("cumulative", "adjacent")) {
          geTerms <- c(geTerms, rval[!subsGe & !subsCe])
      } else if (family$family %in% c("gaussian")) {
          if (any(subsCe) | any(subsGe))
              warning("'ce()' and 'ge()' terms ",
                      "are ignored if family = '",
                      family$family, "'")
          ceTerms <- c(ceTerms, geTerms, rval[!subsGe & !subsCe])
          subsCe <- rep(TRUE, length(ceTerms))
          subsGe <- rep(FALSE, length(rval))
          geTerms <- NULL          
      } else {
          ceTerms <- c(ceTerms, rval[!subsGe & !subsCe])
      }
      ## old
      ## geTerms <- c(geTerms, rval[!subsGe & !subsCe]) # ? better
      
      rval <- list(paste(ceTerms, collapse = "+"),
                   paste(geTerms, collapse = "+"))

    } else {

      rval <- list("", "")
      
    }
    
    ## set intercepts
    if (family$family %in% c("baseline", "gaussian"))
        x$intercept <- "ce"    
    if (x$intercept == "none") {
      int <- list(switch(x$type, fe = "-1", re = "-1", vc = "1"),
                  switch(x$type, fe = "1", re = "-1", vc = "1"))
    } else if (x$intercept == "ge") {
      ## notice that this option is disabled for 'fe' and 'vc' terms
      int <- list(switch(x$type, fe = "1", re = "-1", vc = "1"),
                  switch(x$type, fe = "1", re = "1",
                         vc = paste("Node", x$name, sep = "")))
    } else if (x$intercept == "ce") {
      int <- list(switch(x$type, fe = "1", re = "1",
                         vc = paste("Node", x$name, sep = "")),
                  switch(x$type, fe = "1", re = "-1", vc = "1"))
    }
    
    rval <- lapply(1L:2L, function(i) {
      if (rval[[i]] == "") return(int[[i]])
      if (int[[i]] != "1") return(paste(int[[i]], rval[[i]], sep = "+"))
      return(rval[[i]])
    })
    
    if (x$type == "re" && family$family == "cumulative" && x$intercept == "ce")
      stop("category-specific random effects are not ",
           "available for the cumulative model.")

    if (x$type %in% c("vc", "re") && (rval[[1L]] == "-1" & rval[[2L]] == "-1"))
      stop("the term '", x$type, "' is misspecified any may be dropped from ",
           "the specification")
    
    ## set formulas
    rval <- lapply(paste("~", rval, sep = "") , as.formula)
    
    if (x$type == "re") {
      int <- sapply(rval, function(x) attr(terms(x), "intercept"))
      if (all(int == 1L)) rval[[1L]] <- update(rval[[1L]], ~ . - 1)
    }
    
    ## incorporate 'Node'
    if (x$type == "vc") {

      ## add 'Node' to all terms
      addNode <- function(x, name) {
        terms <- attr(terms(x, keep.order = TRUE), "term.labels")
        termsNew <- paste("Node", name, ":", terms, sep = "")
        if (length(terms) > 0L) {
          x <- update(x, paste("~.-",paste(terms, collapse = "-")))
          x <- update(x, paste("~.+",paste(termsNew, collapse = "+")))
        }
        return(x)
      }
      rval <- lapply(rval, addNode, name = x$name)
    }

  } else {

    if (x$type == "vc") {
      terms <- attr(terms(x$eta, keep.order = TRUE), "term.labels")
      if (length(terms) > 0L) terms <- paste("Node", x$name, ":", terms, sep = "")
      intercept <- x$intercept
      if (x$intercept != "none")
        terms <- c(paste("Node", x$name, sep = ""), terms)
      if (length(terms) == 0L) terms <- "1"
      rval <- list(~1, as.formula(paste("~", paste(terms, collapse = "+"))))
    } else {
      rval <- list(~1, x$eta)
    }
  }
  environment(rval[[1L]]) <- env
  environment(rval[[2L]]) <- env
  names(rval) <- c("ce", "ge")
  return(rval)
}


## --------------------------------------------------------- #
##' Utility function for \code{\link{vcrpart_formula}}
##' 
##' Constructs a formula of conditioning variables, which are
##' the moderators for \code{\link{vc}} terms and the grouping
##' factors for \code{\link{re}} terms.
##'
##' @param x       a formula.
##' @param family  the model family, e.g., cumulative()
##' @param env     environment for evaluating the formula.
##'
##' @return A formula.
##'
##' @details Used in \code{vcrpart_formula}.
## --------------------------------------------------------- #

vcrpart_formula_cond <- function (x, family, env) {
  rval <- x$cond
  environment(rval) <- env
  return(rval)
}


## --------------------------------------------------------- #
##' Construct a formula for \code{\link{tvcm}}, \code{\link{fvcm}}
##' and \code{\link{olmm}} calls
##' 
##' Evaluates input formulas for \code{\link{tvcm}},
##' \code{\link{fvcm}} and \code{\link{olmm}} and returns a
##' list with component wise linear formulas for constructing
##' model matrices.
##'
##' @param formula the original formula.
##' @param family  the model family, e.g., cumulative()
##' @param env     environment for evaluating the formula.
##'
##' @return A list with sublists of formulas.
##'
##' @details Used in \code{\link{fvcm}}, \code{\link{predict.fvcm}},
##' \code{\link{olmm}}, \code{\link{predict.olmm}},
##' \code{\link{tvcm}}, \code{tvcm_get_node},
##' \code{\link{prune.tvcm}}.
## --------------------------------------------------------- #

vcrpart_formula <- function(formula, family = cumulative(),
                            env = parent.frame()) {

  types <- c("fe", "vc", "re")
  terms <- terms(formula, specials = types, keep.order = TRUE)
  yName <- deparse((formula)[[2L]])
  termLabs <- attr(terms, "term.labels")
  termFac <- attr(terms, "factors")
  type <- rep.int("NA", length(termLabs))
  if (length(type) > 0L)
    for (tp in types)
      type[colSums(termFac[attr(terms, "specials")[[tp]],,drop=FALSE]) > 0] <- tp 
  
  getInt <- function(x) eval(parse(text = x))$intercept
  getNPred <- function(x) length(all.vars(eval(parse(text = x))$eta))
  getTerms1 <- function(x, which) {
    x <- eval(parse(text = x))[[which]]
    attr(terms(formula(paste("~", x)), keep.order = TRUE), "term.labels")
  }
  
  ## 'fe' terms
  if (any(type == "fe")) {
    feInt <- unique(unlist(lapply(termLabs[type == "fe"], getInt)))
  } else {
    feInt <- "ce"
  }
  if (length(feInt) > 1L) feInt <- "none"
  feInt <- feInt == "ce"
  feInt <- feInt & attr(terms(formula), "intercept") == 1L
  termLabs[type == "NA"] <- paste("fe(", termLabs[type == "NA"], ")", sep = "")
  type[type == "NA"] <- "fe"  
  feTerms <- unlist(lapply(termLabs[type == "fe"], getTerms1, which = "eta"))
  if (length(feTerms) == 0) feTerms <- "1"
  feTerms <- paste("fe(", paste(feTerms, collapse = " + "),
                   ", intercept = ", feInt, ")", sep = "")
  feTerms <- eval(parse(text = feTerms))
  feTerms$cond <- vcrpart_formula_cond(feTerms, family, env)
  feTerms$eta <- vcrpart_formula_eta(feTerms, family, env)
  
  ## 'vc' terms
  if (any(type == "vc")) {
    vcTerms <- termLabs[type == "vc"]
    vcInt <- unlist(lapply(vcTerms, getInt))
    if (inherits(family, "family")) vcInt[vcInt == "ge"] <- "ce" 
    nPred <- unlist(lapply(vcTerms, getNPred))
    if (feTerms$intercept == "none") {
      direct <- vcInt == "ce"
      direct <- direct &
        unlist(lapply(vcTerms, function(fTerm) !"(Intercept)" %in%
                      eval(parse(text = fTerm))$nuisance))
      subs <- which(direct)[1L]
    } else {
      subs <- c()
    }
    direct <- rep.int(FALSE, length(vcTerms)); direct[subs] <- TRUE;
    
    if (any(direct)) {
      ord <- order(as.integer(!direct))
      direct <- direct[ord]
      vcTerms <- vcTerms[ord]
    }
    vcTerms <- lapply(vcTerms, function(fTerm) eval(parse(text = fTerm)))
    names(vcTerms) <- LETTERS[1:length(vcTerms)]
    for (pid in seq_along(vcTerms)) {
      vcTerms[[pid]]$name <- names(vcTerms)[pid]
      vcTerms[[pid]]$cond <- vcrpart_formula_cond(vcTerms[[pid]], family, env)
      vcTerms[[pid]]$eta <-
        vcrpart_formula_eta(vcTerms[[pid]], family, env)
      subsInt <- vcTerms[[pid]]$nuisance == "(Intercept)"
      vcTerms[[pid]]$nuisance[subsInt] <- paste("Node", LETTERS[pid], sep = "")
      vcTerms[[pid]]$nuisance[!subsInt] <-
        paste("Node", LETTERS[pid], ":", vcTerms[[pid]]$nuisance[!subsInt], sep = "")
      vcTerms[[pid]]$direct <- direct[pid]
    }
  } else {
    vcTerms <- NULL
  }
  
  ## 're' terms
  if (sum(type == "re") > 1L) {
    reCond <- unique(unlist(lapply(termLabs[type == "re"], getTerms1, which = "cond")))
    if (length(reCond) > 1L) stop("non-unique subject factor in 're' terms.")
    reInt <- unique(unlist(lapply(termLabs[type == "re"], getInt)))
    if (length(reInt) == 0L) feInt <- "ge"
    if (length(reInt) > 1L) stop("non-unique 'intercept' in 're' terms.")
    reTerms <- unlist(lapply(termLabs[type == "re"], getTerms1, which = "eta"))
    reTerms <- paste("re(", paste(reTerms, collapse = " + "), "|",
                     reCond, ", intercept = '", reInt, "')", sep = "")
  } else {
    reTerms <- termLabs[type == "re"]
    if (length(reTerms) == 0L) reTerms <- NULL
  }
  if (!is.null(reTerms)) {
    reTerms <- eval(parse(text = reTerms))
    reTerms$cond <- vcrpart_formula_cond(reTerms, family, env)
    reTerms$eta <- vcrpart_formula_eta(reTerms, family, env)
  }

  rval <- list(fe = feTerms,
               vc = vcTerms,
               re = reTerms)
  
  ## add original formula
  rval$original <- as.formula(formula, env = env)

  ## formula with all predictors
  getTerms2 <- function(x) {
    FUN <- function(x) {
      rval <- c(rownames(attr(terms(x$eta$ce), "factors")),
                rownames(attr(terms(x$eta$ge), "factors")),
                rownames(attr(terms(x$cond), "factors")))
      rval <- gsub("Node[A-Z]+:", "", rval)
      rval <- rval[rval != paste("Node", x$name, sep = "")]
    }
    if ("type" %in% names(x)) {
      rval <- rval <- FUN(x)
    } else {
      rval <- unlist(lapply(x, FUN))
    }
    return(rval)
  }
  allTerms <- unique(unlist(lapply(rval[1L:3L], getTerms2)))
  if (length(allTerms) == 0L) allTerms <- "1"
  fAll <- paste(yName, paste(allTerms, collapse = "+"), sep = "~")
  fAll <- as.formula(fAll, env = env)
  rval$all <- fAll
  return(rval)
}


## --------------------------------------------------------- #
##' Delete the environments of formulas in
##' \code{\link{vcrpart_formula}}
##'
##' @param formList a list of formulas as produced by
##'    \code{\link{vcrpart_formula}}.
##'
##' @return A list of formulas.
## --------------------------------------------------------- #

vcrpart_formula_delEnv <- function(formList) {
  if (!is.null(formList$fe)) {
    environment(formList$fe$eta$ce) <- NULL
    environment(formList$fe$eta$ge) <- NULL
  }
  for (pid in seq_along(formList$vc)) {
    environment(formList$vc[[pid]]$eta$ce) <- NULL
    environment(formList$vc[[pid]]$eta$ge) <- NULL
    environment(formList$vc[[pid]]$cond$ce) <- NULL
    environment(formList$vc[[pid]]$cond$ge) <- NULL
  }
  if (!is.null(formList$fe)) {
    environment(formList$fe$eta$ce) <- NULL
    environment(formList$fe$eta$ge) <- NULL
  }
  environment(formList$original) <- NULL
  environment(formList$all) <- NULL
  return(formList)
}
