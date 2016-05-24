### from stats for drop1
safe_pchisq <- function(q, df, ...){
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}

extractAIC.DirichletRegModel <- function(fit, scale = 0, k = 2, ...){
  n <- nobs(fit)
  npar <- fit$npar
  dev <- -2*fit$logLik
  c(npar, dev + k * npar)
}



drop1.DirichletRegModel <- function(
  object,
  scope,
  test = c("LRT", "none"),
  k = 2,
  ...
){

### print caveat
  if(!exists("._DirichletReg_drop1_warning", where = ".GlobalEnv")){
    .GlobalEnv$._DirichletReg_drop1_warning <- TRUE

    if(interactive()) writeLines(paste(rep("-", getOption("width")), collapse = ""))

    writeLines(strwrap(paste(
      "CAVEAT: drop1() is still an experimental feature.",
      "If you plan to use this function, please double-check results, e.g., by comparing two models using anova()."
    , collapse = ""), width = getOption("width"), exdent = 8L))

    writeLines(strwrap("Please don't hesitate to report any problems/doubts: marco.maier@wu.ac.at", width = getOption("width"), exdent = 8L, indent = 8L))

    if(interactive()) writeLines(paste(rep("-", getOption("width")), collapse = ""))
  }
###

  if(!missing(scope)) stop("scope not implemented yet.")

  if(!missing(test) && test == "Chisq"){
    test <- "LRT"
    message('use test = "LRT" instead of "Chisq"')
    test <- match.arg(test, c("LRT", "none"))
  } else if(missing(test)){
    test <- "LRT"
  }

  x <- model.matrix(object)
  n <- nobs(object)

  Formula <- as.Formula(formula(object))

  if(object$parametrization == "common"){
    asgn <- lapply(x[["X"]], attr, which = "assign")
    tl <- lapply(seq_len(length(Formula)[2L]), function(f_index){ attr(terms(Formula, rhs = f_index), "term.labels") })
  } else {
    asgn <- list(
      attr(object$X[[1L]], which = "assign"),
      attr(x[["Z"]], which = "assign")
    )
    tl <- list(
      attr(terms(Formula, rhs = 1L), "term.labels"),
      attr(terms(Formula, rhs = 2L), "term.labels")
    )
  }

  # define scope for variables to drop
  if(missing(scope)){
    scope <- lapply(seq_along(tl), function(f_index){
      drop.scope(formula(Formula, rhs = f_index))
    })
  } else {
stop("not implemented yet!")
    if(!is.character(scope)){
      scope <- attr(terms(update.formula(object, scope)),  "term.labels")
    }
    if(!all(match(scope, tl, 0L) > 0L)){
      stop("scope is not a subset of term labels")
    }
  }

  ndrop <- lapply(seq_along(scope), function(i){ match(scope[[i]], tl[[i]]) })

###
  # naming stuff
  if(object$parametrization == "common"){
    names(asgn) <- names(tl) <- names(scope) <- names(ndrop) <- attr(object$Y, "dim.names")
  } else {
    names(asgn) <- names(tl) <- names(scope) <- names(ndrop) <- c("Mean", "Precision")
  }
###

  ns <- lapply(scope, length)
  chisq <- -2*object$logLik
  dfs <- lapply(ns, numeric)
  dev <- lapply(ns, numeric)

#  y <- object$Y
#  wt <- object$weights

### refit and save values
  for(comp in seq_along(ns)){
    for(subterm in seq_len(ns[[comp]])){
#cat("\ntrying ...", comp, (tl[[comp]])[(ndrop[[comp]][subterm])]); flush.console()
      z <- update(object,
        as.formula(paste0(".~", paste(rep(".|", comp - 1L), collapse=""), ".-", (tl[[comp]])[(ndrop[[comp]][subterm])]))
      )
      dfs[[comp]][subterm] <- z$npar
      dev[[comp]][subterm] <- -2*z$logLik
#cat(" - DONE!\n")
    }
  }

### add components in front of terms
  if(object$parametrization == "alternative"){
    if(length(scope[[1L]]) == 0L){ scope[1L] <- list(NULL) } else { scope[[1L]] <- paste("Mean:", scope[[1L]]) }
    if(length(scope[[2L]]) == 0L){ scope[2L] <- list(NULL) } else { scope[[2L]] <- paste("Prec:", scope[[2L]]) }
  } else if(object$parametrization == "common"){
    for(i in seq_along(scope)){
      if(length(scope[[i]]) == 0L){ scope[i] <- list(NULL) } else { scope[[i]] <- paste0(names(asgn)[i], ": ", scope[[i]]) }
    }
  } else { stop("unexpected error!") }

  scope <- c("<none>",    unlist(scope))
  dfs   <- c(object$npar, unlist(dfs))
  dev   <- c(chisq,       unlist(dev))
  aic <- dev + k * dfs
  dfs <- dfs[1L] - dfs
  dfs[1L] <- NA

  aic <- aic + (AIC(object, k = k) - aic[1L])

  aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = scope, check.names = FALSE)
  if(all(is.na(aic))) aod <- aod[, -3]

  if(test == "LRT"){
      dev <- pmax(0, dev - dev[1L])
      dev[1L] <- NA
      nas <- !is.na(dev)
      LRT <- "LRT"
      aod[, LRT] <- dev
      dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
      aod[, "Pr(>Chi)"] <- dev
  }
  heading <- c("Single term deletions", "\nModel:", deparse(formula(object)))

  class(aod) <- c("anova", "data.frame")

  attr(aod, "heading") <- heading

  return(aod)

}

#drop1(mm)
