anova.mclogit <- function (object, ..., dispersion = NULL, test = NULL) 
{
  dotargs <- list(...)
  named <- if (is.null(names(dotargs))) 
    rep_len(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named)) 
    warning("the following arguments to 'anova.mclogit' are invalid and dropped: ", 
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.mclogit <- vapply(dotargs, function(x) inherits(x, "mclogit") ,
                       #&!inherits(x,"mclogitRandeff"), 
                   NA)
  dotargs <- dotargs[is.mclogit]
  if (length(dotargs)) 
    return(anova.mclogitlist(c(list(object), dotargs), dispersion = dispersion, 
                         test = test))
  stop("'anova.mclogit' can only be used to compare fitted models")
}

anova.mclogitlist <- function (object, ..., dispersion = NULL, test = NULL) 
{
  responses <- as.character(lapply(object, function(x) {
    deparse(formula(x)[[2L]])
  }))
  sameresp <- responses == responses[1L]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning(gettextf("models with response %s removed because response differs from model 1", 
                     sQuote(deparse(responses[!sameresp]))), domain = NA)
  }
  ns <- sapply(object, function(x) x$N)
  if (any(ns != ns[1L])) 
    stop("models were not all fitted to the same size of dataset")
  nmodels <- length(object)
  if (nmodels == 1) stop("'anova.mclogit' can only be used to compare fitted models")
    
  resdf <- as.numeric(lapply(object, function(x) x$residual.df))
  resdev <- as.numeric(lapply(object, function(x) x$deviance))
  table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, 
                                                            -diff(resdev)))
  variables <- lapply(object, function(x) paste(deparse(formula(x)), 
                                                collapse = "\n"))
  dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", 
                                        "Df", "Deviance"))
  title <- "Analysis of Deviance Table\n"
  topnote <- paste("Model ", format(1L:nmodels), ": ", variables, 
                   sep = "", collapse = "\n")
  if (!is.null(test)) {
    bigmodel <- object[[order(resdf)[1L]]]
    df.dispersion <- Inf
    table <- stat.anova(table = table, test = test, scale = 1, 
                        df.scale = df.dispersion, n = bigmodel$N)
  }
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                          "data.frame"))
}


