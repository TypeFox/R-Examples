is.character.or.factor <- function(x) {
  is.character(x) | is.factor(x)
}

is.numeric.and.not.surv <- function(x) {
  is.numeric(x) & !is.Surv(x)
}

##' Cross x and y
##'
##' @param x x
##' @param y y
##' @param funs funs
##' @param ... \dots
##' @param cum cum
##' @param margin margin
##' @param addmargins addmargins
##' @param useNA useNA
##' @param propNA propNA
##' @param revert whether to regroup factors or numeric variables when crossing factor with numeric variables
##' @param method method
##' @param times times
##' @param followup followup
##' @param test test
##' @param test.tabular test.tabular
##' @param test.summarize test.summarize
##' @param test.survival 
##' @param show.test show.test
##' @param plim plim
##' @param show.method show.method
##' @param label display label? (using \code{Hmisc:::label})
##' @author David Hajage
##' @keywords internal
cross <- function(x, y = NULL, funs = c(mean, sd, quantile, n, na), ..., cum = FALSE, margin = 0:2, addmargins = FALSE, useNA = c("no", "ifany", "always"), propNA = TRUE, revert = FALSE, method = c("pearson", "kendall", "spearman"), times = NULL, followup = FALSE, test = FALSE, test.tabular = test.tabular.auto, test.summarize = test.summarize.auto, test.survival = test.survival.logrank, show.test = display.test, plim = 4, show.method = TRUE, label = FALSE) {
  if (!is.character(funs)) {
    funs <- as.character(as.list(substitute(funs)))
    funs <- funs[funs != "c" & funs != "list"]
  }

  results <- "What?"
  
  if (!is.null(x) & !is.null(y)) {
    if (all(sapply(x, is.numeric.and.not.surv)) & all(sapply(y, is.character.or.factor))) {
      results <- summarize.data.frame.by(x, y, funs = funs, ..., addmargins = addmargins, useNA = useNA, revert = revert, test = test, test.summarize = test.summarize, show.test = show.test, plim = plim, show.method = show.method, label = label)
    }
    if (all(sapply(y, is.numeric.and.not.surv)) & all(sapply(x, is.character.or.factor))) {
      results <- summarize.data.frame.by(y, x, funs = funs, ..., addmargins = addmargins, useNA = useNA, revert = revert, test = test, test.summarize = test.summarize, show.test = show.test, plim = plim, show.method = show.method, label = label)
    }
    if (all(sapply(x, is.Surv)) & all(sapply(y, is.character.or.factor))) {
      results <- survival.data.frame.by(x, y, times = times, followup = followup, test = test, test.survival = test.survival, show.test = show.test, plim = plim, show.method = show.method, label = label)
    }
    if (all(sapply(y, is.Surv)) & all(sapply(x, is.character.or.factor))) {
      results <- survival.data.frame.by(y, x, times = times, followup = followup, test = test, test.survival = test.survival, show.test = show.test, plim = plim, show.method = show.method, label = label)
    }
    if (all(sapply(x, is.character.or.factor)) & all(sapply(y, is.character.or.factor))) {
      results <- tabular.data.frame(x, y, margin = margin, addmargins = addmargins, useNA = useNA, propNA = propNA, test = test, test.tabular = test.tabular, show.test = show.test, plim = plim, show.method = show.method, label = label)
    }
    if (all(sapply(x, is.numeric.and.not.surv)) & all(sapply(y, is.numeric.and.not.surv))) {
      results <- correlation.data.frame(x, y, method = method)
    }
  } else if (is.null(y)) {
    if (all(sapply(x, is.character.or.factor))) {
      results <- freq.data.frame(x, addmargins = addmargins, useNA = useNA, propNA = propNA, cum = cum, label = label)
    }
    if (all(sapply(x, is.numeric.and.not.surv))) {
      results <- summarize.data.frame(x, funs = funs, ..., label = label)
    }
    if (all(sapply(x, is.Surv))) {
      results <- survival.data.frame(x, times = times, followup = followup, label = label)
    }
  } else if (is.null(x)) {
    if (all(sapply(y, is.character.or.factor))) {
      results <- freq.data.frame(y, addmargins = addmargins, useNA = useNA, propNA = propNA, cum = cum, label = label)
    } 
    if (all(sapply(y, is.numeric.and.not.surv))) {
      results <- summarize.data.frame(y, funs = funs, ..., label = label)
    }
    if (all(sapply(y, is.Surv))) {
      results <- survival.data.frame(y, times = times, followup = followup, label = label)
    }
  }

  attr(results, "split_type") <- NULL
  attr(results, "split_labels") <- NULL
  
  results
}

##' Cross variables in a list
##'
##' @param l 
##' @param funs funs
##' @param ... \dots
##' @param cum cum
##' @param margin margin
##' @param addmargins addmargins
##' @param useNA useNA
##' @param propNA 
##' @param revert whether to regroup factors or numeric variables when crossing factor with numeric variables
##' @param method method
##' @param times times
##' @param followup followup
##' @param test test
##' @param test.summarize test.summarize
##' @param test.survival 
##' @param test.tabular test.tabular
##' @param show.test show.test
##' @param plim plim
##' @param show.method show.method 
##' @param label label
##' @author David Hajage
##' @keywords internal
cross_list <- function(l, funs = c(mean, sd, quantile, n, na), ..., cum = FALSE, margin = 0:2, addmargins = FALSE, useNA = c("no", "ifany", "always"), propNA = TRUE, revert = FALSE, method = c("pearson", "kendall", "spearman"), times = NULL, followup = FALSE, test = FALSE, test.summarize = test.summarize.auto, test.survival = test.survival.logrank, test.tabular = test.tabular.auto, show.test = display.test, plim = 4, show.method = TRUE, label = FALSE) {

  if (!is.character(funs)) {
    funs <- as.character(as.list(substitute(funs)))
    funs <- funs[funs != "c" & funs != "list"]
  }
  
  x <- l[[1]]
  if (length(l) == 2) {
    y <- l[[2]]
  } else {
    y <- NULL
  }

  cross(x = x, y = y, funs = funs, ..., cum = cum, margin = margin, addmargins = addmargins, useNA = useNA, propNA = propNA, revert = revert, method = method, times = times, followup = followup, test = test, test.summarize = test.summarize, test.tabular = test.tabular, show.test = show.test, plim = plim, show.method = show.method, label = label)
}

##' Regroup factors with factors, and numerical variables with numerical variables
##'
##' @param vars vars
##' @param numdata numdata
##' @param catdata catdata
##' @param survdata survdata
##' @author David Hajage
##' @keywords internal
regroup <- function(vars, numdata, catdata, survdata) {
  vars <- lapply(vars, function(x) remove_blank(elements(x)))

  
  results <- unique(unlist(lapply(vars, function(x) {
    numvars <- x[x %in% numdata]
    catvars <- x[x %in% catdata]
    survvars <- x[x %in% survdata]
    dotvars <- x[x == "."]
    xx <- c(if (length(numvars) > 1) paste("cbind(", paste(numvars, collapse = ","), ")", sep = "") else numvars,
            if (length(catvars) > 1) paste("cbind(", paste(catvars, collapse = ","), ")", sep = "") else catvars,
            if (length(survvars) > 1) paste("cbind(", paste(survvars, collapse = ","), ")", sep = "") else survvars,
            if (length(dotvars) >= 1) ".")
    xx[xx != "cbind()"]
  })))
  
  if (length(results) == 0)
    results <- "."
  results
}

##' Remix and describe.
##'
##' A quick and easy function for describing datasets.
##'
##' @param formula a formula (see Details).
##' @param data a data.frame.
##' @param funs functions for describing numeric variable.   Can be
##' \code{c(fun1, fun2, fun3)} or   \code{c("fun1", "fun2", "fun3")}
##' or a list.
##' @param ... further arguments (all passed to funs), for example
##'{na.rm = TRUE}
##' @param cum should cumulated frequencies be reported?
##' @param margin index, or vector of indices to generate proportion
##' in frequency tables (0: cell, 1: row, 2: col).
##' @param addmargins whether to add margins
##' @param useNA whether to include NA as a level (factor)
##' @param propNA whether to include NA in proportion calculation
##' @param revert whether to regroup factors or numeric variables when crossing factor with numeric variables
##' @param method a character string indicating which correlation
##' coefficient is to be   used. One of \code{"pearson"},
##' \code{"kendall"}, or \code{"spearman"}, can be abbreviated.
##' @param times vector of times (see \code{?summary.survival}
##' un package \code{survival})
##' @param followup whether to display follow-up time
##' @param test whether to perform tests
##' @param test.summarize a function of two arguments (continuous
##' variable and grouping variable) used to compare continuous
##' variable, that return a list of two components : \code{p.value}
##' and \code{method} (the test name). See \code{test.summarize.auto},
##' \code{test.summarize.kruskal},
##' \code{test.summarize.oneway.equalvar}, or
##' \code{test.summarize.unequalvar} for example of such
##' functions. Users can provide their own function.
##' @param test.survival a function of one argument (a formula) used
##' to compare survival estimations, that returns the same components
##' as created by \code{test.summarize}. See
##' \code{test.survival.logrank}. Users can provide their own
##' function.
##' @param test.tabular a function of three arguments (two categorical
##' variables and a logical \code{na}) used to test association
##' between two factors, that returns the same components as created
##' by \code{test.summarize}. See \code{test.tabular.auto} and
##' \code{test.tabular.fisher}. Users can provide their own function.
##' @param show.test a function used to display the test result. See
##' \code{display.test}.
##' @param plim number of digits for the p value
##' @param show.method should display the test name?
##' @param label whether to display labels of variables (using
##' \code{label} in package \code{Hmisc})
##' @note
##'   The formula has the following format: \code{x_1 + x_2 + ... ~ y_1 + y_2 + ...}
##'
##'   There are a couple of special variables: \code{...} represents all
##'   other variables not used in the formula and \code{.} represents no
##'   variable, so you can do \code{formula = var1 ~ .}.
##'
##'   If \code{var1} is numeric, \code{var1 ~ .} produce a summary
##'   table using \code{funs}. If \code{var1} is a factor, \code{var1 ~
##'   .} produce a frequency table. If \code{var1} is of class
##'   \code{Surv}, \code{var1 ~ .} produce a table with the estimates of
##'   survival at \code{times}. If \code{var1} is numeric and
##'   \code{var2} is numeric, \code{var1 ~ var2} gives correlation. if
##'   \code{var1} is numeric and \code{var2} is a factor, \code{var1 ~
##'   var2} produce a summary table using \code{funs} according to the
##'   levels of \code{var2}. If \code{var1} is a factor and \code{var2}
##'   is a factor, \code{var1 ~ var2} produce a contingency table. If
##'   \code{var1} is of class \code{Surv} and \code{var2} is a factor,
##'   \code{var1 ~ var2} produce a table with the estimates of survival
##'   for each level of \code{var2}.
##'
##'   You can group several variables of the same type (numeric or factor)
##'   together with \code{cbind(var1, var2, var3)}, they will be grouped in the
##'   same table. \code{cbind(...)} works (ie regroups all variables of the same
##'   type).
##'
##' 
##' @return
##'   A remix object, basically a list with descriptive tables. It uses
##'   \code{ascii} package for printing output, and can be use with
##'   \code{ascii} function.
##' @author David Hajage, inspired by the design and the code of
##'   \code{summary.formula} (\code{Hmisc} package, FE Harrell) and
##'   \code{cast} (\code{reshape} package, H Wickham).
##' @seealso \code{cast} (reshape) and \code{summary.formula} (Hmisc).
##' @examples
##' parwidth <- getOption("width")
##' options(width = 100)
##'
##' library(remix)
##' remix(data = iris)
##' remix(cbind(...) ~ ., iris[, sapply(iris, is.numeric)], funs = c(median, mad, min, max))
##' remix(cbind(Sepal.Length, I(Sepal.Width^2)) ~ Species, iris, funs = quantile, probs = c(1/3, 2/3))
##' remix(Sepal.Length + Sepal.Width ~ Petal.Length + Petal.Width, iris)
##' remix(cbind(Sepal.Length, Sepal.Width) ~ cbind(Petal.Length, Petal.Width), iris)
##' remix(... ~ ., esoph, cum = TRUE)
##' remix(alcgp ~ tobgp, esoph, cum = TRUE)
##' remix(Surv(time, status) ~ x, data = aml, times = seq(0, 120, 12))
##' 
##' options(width = parwidth)
##' @keywords univar
##' @export
remix <- function(formula = cbind(...) ~ ., data = NULL, funs = c(mean, sd, quantile, n, na), ..., cum = FALSE, margin = 0:2, addmargins = FALSE, useNA = c("no", "ifany", "always"), propNA = TRUE, revert = FALSE, method = c("pearson", "kendall", "spearman"), times = NULL, followup = FALSE, test = FALSE, test.summarize = test.summarize.auto, test.survival = test.survival.logrank, test.tabular = test.tabular.auto, show.test = display.test, plim = 4, show.method = TRUE, label = FALSE) {
  
  if (is.formula(formula))
    formula <- deparse(formula, 500)
  
  if (!is.character(funs)) {
    funs <- as.character(as.list(substitute(funs)))
    funs <- funs[funs != "c" & funs != "list"]
  }

  varnames <- names(data)
  parsed <- parse_formula(formula, varnames)
  
  data <-   parse_data(expand_formula(formula, varnames), data)
  names(data) <- remove_blank(names(data))
  varform <- names(data)

  numdata <- varform[sapply(data, function(x) is.numeric(x) & !is.Surv(x))]
  catdata <- varform[sapply(data, is.character.or.factor)]
  survdata <- varform[sapply(data, is.Surv)]
  
  parsed$left <- regroup(parsed$left, numdata, catdata, survdata)
  parsed$right <- regroup(parsed$right, numdata, catdata, survdata)
  
  eg <- expand.grid(parsed$left, parsed$right)
  
  ## if (all(parsed$by == ".")) {
    comb <- lapply(apply(eg, 1, list), function(x) {
      y <- unlist(x)
      y <- y[y != "."]
      ## y <- sub("(cbind\\()(.*)(\\))", "\\2", y)
      ## lapply(y, function(z) data[, strsplit(z, ",")[[1]], drop = FALSE])})
      lapply(y, function(z) data[, remove_blank(elements(z)), drop = FALSE])})
    
    results <- llply(comb, cross_list, funs = funs, ..., cum = cum, margin = margin, addmargins = addmargins, useNA = useNA, propNA = propNA, revert = revert, method = method, times = times, followup = followup, test = test, test.summarize = test.summarize, test.tabular = test.tabular, show.test = show.test, plim = plim, show.method = show.method, label = label)
    names(results) <- apply(eg, 1, paste, collapse = " ~ ")
  ## }
  
  ## if (any(parsed$by != ".")) {
  ##   databy <- list()
  ##   for (i in 1:length(parsed$by)) {
  ##     ldata <- dlply(data, parsed$by[i], function(x){x})
  ##     attr(ldata, "split_type") <- NULL
  ##     attr(ldata, "split_labels") <- NULL
      
  ##     databy <- c(databy, list(ldata))
  ##   }
  ##   names(databy) <- parsed$by

  ##   combby <- lapply(databy, function(o) lapply(o, function(p) lapply(apply(eg, 1, list), function(x) {
  ##     y <- unlist(x)
  ##     y <- y[y != "."]
  ##     y <- sub("(cbind *\\()(.*)(\\))", "\\2", y)
  ##     lapply(y, function(z) p[, strsplit(z, ",")[[1]], drop = FALSE])})))

  ##   results <- lapply(combby, function(x) lapply(x, function(comb) {
  ##     res <- llply(comb, cross_list, funs = funs, cum = cum, margin = margin, addmargins = addmargins, useNA = useNA, propNA = propNA, revert = revert, method = method)
  ##     names(res) <- apply(eg, 1, paste, collapse = " ~ ")
  ##     res
  ##   }))
  ##   if (length(results) == 1)
  ##     results <- results[[1]]
  ## }
  
  class(results) <- c("remix")
  attr(results, "formula") <- formula
  attr(results, "left") <- parsed$left
  attr(results, "right") <- parsed$right
  ## attr(results, "by") <- parsed$by
  
  attr(results, "data") <- data
  results
}

##' Ascii for remix object.
##'
##' Ascii method for remix object.
##'
##' @export
##' @method ascii remix
##' @import ascii
##' @param x a remix object
##' @param caption.level see \code{?ascii} in \code{ascii} package
##' @param format see \code{?ascii} in \code{ascii} package
##' @param digits see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii} (all except \code{caption}
##'    which has no effect)
##' @author David Hajage
##' @keywords univar
ascii.remix <- function(x, caption.level = c("s", "e", "m"), format = "nice", digits = 2, ...) {
  caption.level <- rep(caption.level, length = 3)
  caption.level1 <- caption.level[1]
  caption.level2 <- caption.level[2]
  caption.level3 <- caption.level[3]
  
  xx <- list()
  ## if (all(attr(x, "by") == ".")) {
    captions <- names(x)
    for (i in 1:length(x)) {
      xx[[i]] <- ascii(x[[i]], caption = captions[i], caption.level = caption.level1, format = format, digits = digits, ...)
    }
  ## } else if (length(attr(x, "by")) == 1) {
  ##   captions1 <- names(x)
  ##   captions2 <- names(x[[1]])
  ##   for (i in 1:length(x)) {
  ##     asc.cap1 <- ascii(list(NULL), caption = captions1[i], caption.level = caption.level1)
  ##     xx[[paste("obj", i, sep = "")]] <- asc.cap1
  ##     for (j in 1:length(x[[i]])) {
  ##       xx[[paste("obj", i, j, sep = "")]] <- ascii(x[[i]][[j]], caption = captions2[j], caption.level = caption.level2, format = format, digits = digits, ...)
  ##     }
  ##   }
  ## } else if (length(attr(x, "by")) > 1) {
  ##   captions1 <- names(x)
  ##   captions3 <- names(x[[1]][[1]])
  ##   for (i in 1:length(x)) {
  ##     asc.cap1 <- ascii(list(NULL), caption = captions1[i], caption.level = caption.level1)
  ##     xx[[paste("obj", i, sep = "")]] <- asc.cap1
  ##     for (j in 1:length(x[[i]])) {
  ##       captions2 <- names(x[[i]])
  ##       asc.cap2 <- ascii(list(NULL), caption = captions2[j], caption.level = caption.level2)
  ##       xx[[paste("obj", i, j, sep = "")]] <- asc.cap2
  ##       for (k in 1:length(x[[i]][[j]])) {
  ##         xx[[paste("obj", i, j, k, sep = "")]] <- ascii(x[[i]][[j]][[k]], caption = captions3[k], caption.level = caption.level3, format = format, digits = digits)
  ##       }
  ##     }
  ##   }
  ## }
  asciiMixed$new(args = xx)
}

##' Print a remix object
##'
##' Print remix object using ascii package
##'
##' @export
##' @import ascii
##' @method print remix
##' @param x a remix object
##' @param type type of output. See \code{?ascii} in \code{ascii} package
##' @param caption.level see \code{?ascii} in \code{ascii} package
##' @param lstyle see \code{?ascii} in \code{ascii} package
##' @param tstyle see \code{?ascii} in \code{ascii} package
##' @param ... other arguments passed to \code{ascii} (all except \code{caption}
##'    which has no effect)
##' @author David Hajage
##' @keywords univar
print.remix <- function(x, type = "rest", caption.level = 1:3, lstyle = "", tstyle = "", ...) {
  print(ascii.remix(x, caption.level = caption.level, lstyle = lstyle, tstyle = tstyle, ...), type = type)
  ## invisible(x)
}

##' Demix
##'
##' Transfrom a remix object into a (list of) data.frame(s).
##'
##' @export
##' @param x a remix object
##'
##' @return
##'   A list of data.frame.
##' @author David Hajage
##' @seealso \code{remix}
##' @examples
##'   x <- remix(... ~ ., esoph, cum = TRUE)
##'   demix(x)
demix <- function(x) {
  ## if (all(attr(x, "by") == ".")) {
    result <- lapply(x, as.data.frame)
  ## } else if (length(attr(x, "by")) == 1) {
  ##   result <- lapply(x, function(x) lapply(x, as.data.frame))
  ## } else if (length(attr(x, "by")) > 1) {
  ##   result <- lapply(x, function(x) lapply(x, function(x) lapply(x, as.data.frame)))
  ## }
  return(result)
}

##' Test if \code{x} is an remix object
##'
##' Test if \code{x} is an remix object
##'
##' @param x a remix object
##' @author David Hajage
##' @keywords internal
is.remix <- function(x)
    inherits(x, "remix")
