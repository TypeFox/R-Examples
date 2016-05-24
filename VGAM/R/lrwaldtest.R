# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












update_default <- function (object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")

  extras <- match.call(expand.dots = FALSE)$...

  if (!missing(formula.)) {
    call$formula <- update_formula(formula(object), formula.)
  }

  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing])
      call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  if (evaluate) eval(call, parent.frame()) else call
}





update_formula <- function (old, new, ...) {


  tmp <-          (update.formula(as.formula(old), as.formula(new)))




  out <- formula(terms.formula(tmp, simplify = TRUE))
  return(out)
}







if (FALSE)
print_anova <- function (x, digits = max(getOption("digits") - 2, 3),
                         signif.stars = getOption("show.signif.stars"),
                         ...) {


  x <- x@Body

  if (!is.null(heading <- attr(x, "heading")))
    cat(heading, sep = "\n")
  nc <- dim(x)[2L]
  if (is.null(cn <- colnames(x)))
    stop("'anova' object must have colnames")
  has.P <- grepl("^(P|Pr)\\(", cn[nc])
  zap.i <- 1L:(if (has.P) nc - 1 else nc)
  i <- which(substr(cn, 2, 7) == " value")
  i <- c(i, which(!is.na(match(cn, c("F", "Cp", "Chisq")))))
  if (length(i))
    zap.i <- zap.i[!(zap.i %in% i)]
  tst.i <- i
  if (length(i <- grep("Df$", cn)))
    zap.i <- zap.i[!(zap.i %in% i)]
  stats::printCoefmat(x, digits = digits, signif.stars = signif.stars, 
      has.Pvalue = has.P, P.values = has.P, cs.ind = NULL, 
      zap.ind = zap.i, tst.ind = tst.i, na.print = "", ...)
  invisible(x)
}



  setGeneric("lrtest", function(object, ...) standardGeneric("lrtest"),
             package = "VGAM")






setClass("VGAManova", representation(
         "Body"        = "data.frame"))



lrtest_vglm <- function(object, ..., name = NULL) {








  cls <- class(object)[1]

  nobs <- function(x) x@misc$nrow.X.vlm


  tlab <- function(x) attr(terms(x), "term.labels")


  if (is.null(name))
    name <- function(x) paste(deparse(formula(x)), collapse = "\n")



  modelUpdate <- function(fm, update) {

    if (is.numeric(update)) {
      if (any(update < 1)) {
        warning("for numeric model specifications all values ",
                "have to be >=1")
        update <- abs(update)[abs(update) > 0]
      }
      if (any(update > length(tlab(fm)))) {
        warning(paste("more terms specified than existent in the model:",
                paste(as.character(update[update > length(tlab(fm))]),
                      collapse = ", ")))
        update <- update[update <= length(tlab(fm))]
      }
      update <- tlab(fm)[update]
    }

    if (is.character(update)) {
      if (!all(update %in% tlab(fm))) {
        warning(paste("terms specified that are not in the model:",
                paste(dQuote(update[!(update %in% tlab(fm))]),
                      collapse = ", ")))
        update <- update[update %in% tlab(fm)]
      }
      if (length(update) < 1) stop("empty model specification")  
      update <- as.formula(paste(". ~ . -",
                           paste(update, collapse = " - ")))
    }

    if (inherits(update, "formula")) {
      update <- update_default(fm, update)
    }
    if (!inherits(update, cls)) {
      warning(paste("original model was of class \"", cls,
                    "\", updated model is of class \"",
                    class(update)[1], "\"", sep = ""))
    }
    return(update)
  }


  objects <- list(object, ...)
  nmodels <- length(objects)
  if (nmodels < 2) {
    objects <- c(objects, . ~ 1)
    nmodels <- 2
  }
  
  no.update <- sapply(objects, function(obj) inherits(obj, cls))
  
  for (i in 2:nmodels) {
    objects[[i]] <- modelUpdate(objects[[i-1]], objects[[i]])
  }


  ns <- sapply(objects, nobs)


  if (any(ns != ns[1])) {
    for (i in 2:nmodels) {
      if (ns[1] != ns[i]) {
        if (no.update[i])
          stop("models were not all fitted to ",
               "the same size of dataset") else {
            commonobs <- row.names(model.frame(objects[[i]])) %in%
                         row.names(model.frame(objects[[i-1]]))
            objects[[i]] <- eval(substitute(update(objects[[i]],
                                 subset = commonobs),
                                 list(commonobs = commonobs)))
            if (nobs(objects[[i]]) != ns[1])
             stop("models could not be fitted to the same size of dataset")
        }
      }
    }
  }

  rval <- matrix(rep(NA_real_, 5 * nmodels), ncol = 5)
  colnames(rval) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(rval) <- 1:nmodels
  
  logLlist <- lapply(objects, logLik)

  dflist <- lapply(objects, df.residual)

  rval[,1] <- unlist(dflist)

  rval[,2] <- unlist(logLlist)

  rval[2:nmodels, 3] <- rval[2:nmodels, 1] - rval[1:(nmodels-1), 1]

  rval[2:nmodels, 4] <- 2 * abs(rval[2:nmodels, 2] - rval[1:(nmodels-1), 2])

  rval[,5] <- pchisq(rval[,4], round(abs(rval[,3])), lower.tail = FALSE)

  variables <- lapply(objects, name)
  title <- "Likelihood ratio test\n"
  topnote <- paste("Model ", format(1:nmodels),
                   ": ", variables, sep = "", collapse = "\n")





  new("VGAManova", Body =
    structure(as.data.frame(rval), heading = c(title, topnote)))
}




setMethod("lrtest", "vglm",
         function(object, ...)
         lrtest_vglm(object = object, ...))












 setMethod("show", "VGAManova",
           function(object)
             getS3method("print", "anova")(object@Body))









use.S3.lrtest <- TRUE
use.S3.lrtest <- FALSE


if (use.S3.lrtest)
lrtest <- function(object, ...) {
  UseMethod("lrtest")
}



if (use.S3.lrtest)
lrtest.formula <- function(object, ..., data = list()) {
  object <- if (length(data) < 1)
    eval(call("lm", formula = as.formula(deparse(substitute(object))),
              environment(object))) else
    eval(call("lm", formula = as.formula(deparse(substitute(object))),
              data = as.name(deparse(substitute(data))),
              environment(data)))
  lrtest.default(object, ...)
}



if (use.S3.lrtest)
lrtest.default <- function(object, ..., name = NULL) {






  cls <- class(object)[1]

  nobs <- function(x) NROW(residuals(x))
  tlab <- function(x) attr(terms(x), "term.labels")
  if (is.null(name))
    name <- function(x) paste(deparse(formula(x)), collapse = "\n")

  modelUpdate <- function(fm, update) {
    if (is.numeric(update)) {
      if (any(update < 1)) {
        warning("for numeric model specifications all values ",
                "have to be >=1")
        update <- abs(update)[abs(update) > 0]
      }
      if (any(update > length(tlab(fm)))) {
        warning(paste("more terms specified than existent in the model:",
                paste(as.character(update[update > length(tlab(fm))]),
                      collapse = ", ")))
        update <- update[update <= length(tlab(fm))]
      }
      update <- tlab(fm)[update]
    }
    if (is.character(update)) {
      if (!all(update %in% tlab(fm))) {
        warning(paste("terms specified that are not in the model:",
                paste(dQuote(update[!(update %in% tlab(fm))]),
                      collapse = ", ")))
        update <- update[update %in% tlab(fm)]
      }
      if (length(update) < 1) stop("empty model specification")  
      update <- as.formula(paste(". ~ . -",
                           paste(update, collapse = " - ")))
    }
    if (inherits(update, "formula")) update <- update(fm, update)
    if (!inherits(update, cls))
      warning(paste("original model was of class \"", cls,
                    "\", updated model is of class \"",
                    class(update)[1], "\"", sep = ""))
    return(update)
  }


  objects <- list(object, ...)
  nmodels <- length(objects)
  if (nmodels < 2) {
    objects <- c(objects, . ~ 1)
 print("objects 1")
 print( objects )
    nmodels <- 2
  }
  
  no.update <- sapply(objects, function(obj) inherits(obj, cls))
 print("no.update")
 print( no.update )
  
  for (i in 2:nmodels)
    objects[[i]] <- modelUpdate(objects[[i-1]], objects[[i]])

 print("objects i")
 print( objects )

  ns <- sapply(objects, nobs)
  if (any(ns != ns[1])) {
    for (i in 2:nmodels) {
      if (ns[1] != ns[i]) {
        if (no.update[i])
          stop("models were not all fitted to ",
               "the same size of dataset") else {
            commonobs <- row.names(model.frame(objects[[i]])) %in%
                         row.names(model.frame(objects[[i-1]]))
 print("commonobs")
 print( commonobs )
            objects[[i]] <- eval(substitute(update(objects[[i]],
                                 subset = commonobs),
                                 list(commonobs = commonobs)))
            if (nobs(objects[[i]]) != ns[1])
             stop("models could not be fitted to the same size of dataset")
        }
      }
    }
  }

  rval <- matrix(rep(NA_real_, 5 * nmodels), ncol = 5)
  colnames(rval) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(rval) <- 1:nmodels
  
  logL <- lapply(objects, logLik)
  rval[,1] <- as.numeric(sapply(logL, function(x) attr(x, "df")))  
  rval[,2] <- sapply(logL, as.numeric)
  rval[2:nmodels, 3] <- rval[2:nmodels, 1] - rval[1:(nmodels-1), 1]
  rval[2:nmodels, 4] <- 2 * abs(rval[2:nmodels, 2] - rval[1:(nmodels-1), 2])
  rval[,5] <- pchisq(rval[,4], round(abs(rval[,3])), lower.tail = FALSE)

  variables <- lapply(objects, name)
  title <- "Likelihood ratio test\n"
  topnote <- paste("Model ", format(1:nmodels),
                   ": ", variables, sep = "", collapse = "\n")

  structure(as.data.frame(rval), heading = c(title, topnote),
            class = c("anova", "data.frame"))
} # End of lrtest.default








if (FALSE)
  setGeneric("waldtest", function(object, ...) standardGeneric("waldtest"),
             package = "VGAM")


if (FALSE)
waldtest <- function(object, ...) {
  UseMethod("waldtest")
}



waldtest_formula <- function(object, ..., data = list()) {

  stop("cannot find waldtest_lm()")

  object <- if (length(data) < 1)
    eval(call("lm", formula = as.formula(deparse(substitute(object))),
         environment(object))) else
    eval(call("lm", formula = as.formula(deparse(substitute(object))),
         data = as.name(deparse(substitute(data))), environment(data)))
 
}




waldtest_default <- function(object, ..., vcov = NULL,
                             test = c("Chisq", "F"), name = NULL) {


  vcov. <- vcov
  cls <- class(object)[1]


  nobs <- function(x) NROW(residuals(x))

  tlab <- function(x) attr(terms(x), "term.labels")

  if (is.null(name))
    name <- function(x) paste(deparse(formula(x)), collapse = "\n")



  modelUpdate <- function(fm, update) {
    if (is.numeric(update)) {
      if (any(update < 1)) {
        warning("for numeric model specifications all values ",
                "have to be >=1")
        update <- abs(update)[abs(update) > 0]
      }
      if (any(update > length(tlab(fm)))) {
        warning(paste("more terms specified than existent in the model:",
                paste(as.character(update[update > length(tlab(fm))]),
                      collapse = ", ")))
        update <- update[update <= length(tlab(fm))]
      }
      update <- tlab(fm)[update]
    }
    if (is.character(update)) {
      if (!all(update %in% tlab(fm))) {
        warning(paste("terms specified that are not in the model:",
                paste(dQuote(update[!(update %in% tlab(fm))]),
                      collapse = ", ")))
        update <- update[update %in% tlab(fm)]
      }
      if (length(update) < 1) stop("empty model specification")  
      update <- as.formula(paste(". ~ . -",
                           paste(update, collapse = " - ")))
    }
    if (inherits(update, "formula")) update <- update(fm, update)
    if (!inherits(update, cls))
      stop(paste("original model was of class \"", cls,
                 "\", updated model is of class \"",
                 class(update)[1], "\"", sep = ""))
    return(update)
  }


  modelCompare <- function(fm, fm.up, vfun = NULL) {
    q <- length(coef(fm)) - length(coef(fm.up))

    if (q > 0) {
      fm0 <- fm.up
      fm1 <- fm
    } else {
      fm0 <- fm
      fm1 <- fm.up
    }
    k <- length(coef(fm1))
    n <- nobs(fm1)

    if (!all(tlab(fm0) %in% tlab(fm1)))
      stop("models are not nested")
    ovar <- which(!(names(coef(fm1)) %in% names(coef(fm0))))

    vc <- if (is.null(vfun)) vcov(fm1) else if (is.function(vfun))
          vfun(fm1) else vfun

    stat <- t(coef(fm1)[ovar]) %*% solve(vc[ovar,ovar]) %*% coef(fm1)[ovar]
    return(c(-q, stat))
  }


  objects <- list(object, ...)
  nmodels <- length(objects)
  if (nmodels < 2) {
    objects <- c(objects, . ~ 1)
    nmodels <- 2
  }
  
  no.update <- sapply(objects, function(obj) inherits(obj, cls))
  
  for (i in 2:nmodels)
    objects[[i]] <- modelUpdate(objects[[i-1]], objects[[i]])

  responses <- as.character(lapply(objects,
               function(x) deparse(terms(x)[[2]])))
  sameresp <- responses == responses[1]
  if (!all(sameresp)) {
    objects <- objects[sameresp]
    warning("models with response ", deparse(responses[!sameresp]),
            " removed because response differs from ", "model 1")
  }

  ns <- sapply(objects, nobs)
  if (any(ns != ns[1])) {
    for (i in 2:nmodels) {
      if (ns[1] != ns[i]) {
        if (no.update[i])
          stop("models were not all fitted to the ",
               "same size of dataset") else {
            commonobs <- row.names(model.frame(objects[[i]])) %in%
                         row.names(model.frame(objects[[i-1]]))
            objects[[i]] <- eval(substitute(update(objects[[i]],
                                 subset = commonobs),
                              list(commonobs = commonobs)))
            if (nobs(objects[[i]]) != ns[1])
             stop("models could not be fitted to the same size of dataset")
          }
      }
    }
  }

  if (nmodels > 2 && !is.null(vcov.) && !is.function(vcov.))
    stop("to compare more than 2 models `vcov.' needs to be a function")

  test <- match.arg(test)
  rval <- matrix(rep(NA_real_, 4 * nmodels), ncol = 4)
  colnames(rval) <- c("Res.Df", "Df", test,
                      paste("Pr(>", test, ")", sep = ""))
  rownames(rval) <- 1:nmodels
  rval[,1] <- as.numeric(sapply(objects, df.residual))
  for (i in 2:nmodels)
    rval[i, 2:3] <- modelCompare(objects[[i-1]], objects[[i]],
                                 vfun = vcov.)
  if (test == "Chisq") {
    rval[,4] <- pchisq(rval[,3], round(abs(rval[,2])), lower.tail = FALSE)
  } else {
    df <- rval[,1]
    for (i in 2:nmodels)
      if (rval[i, 2] < 0)
        df[i] <- rval[i-1, 1]
    rval[, 3] <- rval[, 3] / abs(rval[, 2])
    rval[, 4] <- pf(rval[, 3], abs(rval[, 2]), df, lower.tail = FALSE)
  }


  variables <- lapply(objects, name)
  title <- "Wald test\n"
  topnote <- paste("Model ", format(1:nmodels),
                   ": ", variables, sep = "", collapse = "\n")




  new("VGAManova", Body =
    structure(as.data.frame(rval), heading = c(title, topnote)))
}







if (FALSE)
setMethod("waldtest", "vglm",
         function(object, ...)
         waldtest_vglm(object = object, ...))












