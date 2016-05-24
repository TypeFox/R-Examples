anoas <- function(tab, nd=3, symmetric=FALSE, diagonal=FALSE, ...) {
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(symmetric && nrow(tab) != ncol(tab))
      stop("tab must be a square table for symmetric model")

  if(symmetric && !all(rownames(tab) == colnames(tab)))
      stop("tab must have identical row and column names for symmetric model")

  if(is.na(nd) || nd < 1)
      stop("nd should be at least 1")

  if(symmetric && nd/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!symmetric && nd > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  if(length(dim(tab)) > 2)
      tab <- margin.table(tab, 1:2)


  models <- vector("list", nd + 1)
  names(models) <- c("indep", paste("rc", 1:nd, sep=""))


  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  if(diagonal)
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""



  cat("Fitting independence model...\n")

  args <- list(formula=as.formula(sprintf("Freq ~ %s + %s %s", vars[1], vars[2], diagstr)),
               data=tab, family=poisson)

  models[[1]] <- do.call("gnm", c(args, list(...)))


  cat("Fitting model with 1 dimension...\n")

  models[[2]] <- rc(tab, nd=1, symmetric=symmetric, diagonal=diagonal, ...)

  if(is.null(models[[2]]))
      stop("Model with 1 dimension could not be fitted.")

  if(nd > 1) {
      npar <- if(symmetric) nrow(tab) else nrow(tab) + ncol(tab)

      for(i in 2:nd) {
          if(is.null(models[[i]]))
              stop(sprintf("Model with %i dimensions could not be fitted.", i-1))

          cat(sprintf("Fitting model with %i dimensions...\n", i))

          models[[i+1]] <- rc(tab, nd=i, symmetric=symmetric, diagonal=diagonal,
                              start=c(rep(NA, length(parameters(models[[i]])) - (i-1) * npar),
                                      tail(parameters(models[[i]]), (i-1) * npar), rep(NA, npar)),
                              ...)
      }
  }

  class(models) <- "anoas"
  attr(models, "symmetric") <- symmetric

  models
}

anoasL <- function(tab, nd=3, layer.effect=c("homogeneous.scores", "heterogeneous", "none"),
                   symmetric=FALSE, diagonal=c("none", "heterogeneous", "homogeneous"), ...) {
  layer.effect <- match.arg(layer.effect)
  diagonal <- match.arg(diagonal)

  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

  if(is.na(nd) || nd < 1)
      stop("nd should be at least 1")

  if(symmetric && nd/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!symmetric && nd > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  if(length(dim(tab)) > 3)
      tab <- margin.table(tab, 1:3)


  models <- vector("list", nd + 1)
  names(models) <- c("indep", paste("rcL", 1:nd, sep=""))


  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  if(diagonal == "heterogeneous")
      diagstr <- sprintf("+ %s:Diag(%s, %s) ", vars[3], vars[1], vars[2])
  else if(diagonal == "homogeneous")
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""

  eliminate <- eval(parse(text=sprintf("quote(%s:%s)", vars[1], vars[3])))

  cat("Fitting conditional independence model...\n")

  args <- list(formula=as.formula(sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s",
                                          vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3], diagstr)),
               data=tab, family=poisson, eliminate=eliminate)

  models[[1]] <- do.call("gnm", c(args, list(...)))


  cat("Fitting model with 1 dimension...\n")

  models[[2]] <- rcL(tab, nd=1, layer.effect=layer.effect, symmetric=symmetric, diagonal=diagonal, ...)

  if(is.null(models[[2]]))
      stop("Model with 1 dimension could not be fitted.")

  if(nd > 1) {
      npar <- if(layer.effect == "homogeneous.scores") {
                  if(symmetric) dim(tab)[3] + nrow(tab) else dim(tab)[3] + nrow(tab) + ncol(tab)
              }
              else if (layer.effect == "heterogeneous") {
                  if(symmetric) dim(tab)[3] * nrow(tab) else dim(tab)[3] * (nrow(tab) + ncol(tab))
              }
              else {
                  if(symmetric) nrow(tab) else nrow(tab) + ncol(tab)
              }

      for(i in 2:nd) {
          if(is.null(models[[i]]))
              stop(sprintf("Model with %i dimensions could not be fitted.", i-1))

          cat(sprintf("Fitting model with %i dimensions...\n", i))

          models[[i+1]] <- rcL(tab, nd=i, layer.effect=layer.effect, symmetric=symmetric, diagonal=diagonal,
                               start=c(rep(NA, length(parameters(models[[i]])) - (i-1) * npar),
                                       tail(parameters(models[[i]]), (i-1) * npar), rep(NA, npar)),
                               ...)
      }
  }

  class(models) <- c("anoasL", "anoas")
  attr(models, "symmetric") <- symmetric

  models
}

print.anoas <- function(x, ...) {
    result <- summary.anoas(x)
    print(result, ...)

    invisible(result)
}

summary.anoas <- function(object, ...) {
  df <- sapply(object, function(model) if(!is.null(model)) model$df.residual else NA)
  dev <- sapply(object, function(model) if(!is.null(model)) model$deviance else NA)
  diss <- sapply(object, function(model) {
      if(!is.null(model)) sum(na.omit(abs(c(residuals(model, "response")))))/sum(na.omit(abs(c(fitted(model)))))/2
      else NA })
  bic <- sapply(object, function(model) if(!is.null(model)) model$deviance - log(sum(na.omit(c(model$data)))) * model$df.residual else NA)
  aic <- sapply(object, function(model) if(!is.null(model)) model$deviance - 2 * model$df.residual else NA)

  result <- data.frame(df, dev, dev/dev[1] * 100, diss * 100, bic, aic, c(NA, diff(dev)), c(NA, diff(df)))

  names(result) <- c("Res. Df", "Res. Dev", "Dev./Indep. (%)", "Dissim. (%)", "BIC", "AIC", "Dev.", "Df")

  if(inherits(object, "anoasL"))
      rownames(result) <- c("Conditional indep.", paste("RC(", seq(1, length(object) - 1), ")-L",
                                                        if(attr(object, "symmetric")) " symmetric" else "", sep=""))
  else
      rownames(result) <- c("Indep.", paste("RC(", seq(1, length(object) - 1), ")",
                                            if(attr(object, "symmetric")) " symmetric" else "", sep=""))

  class(result) <- c("summary.anoas", "data.frame")

  result
}

print.summary.anoas <- function(x, digits=1, nsmall=2, scientific=FALSE, ...) {
  print(format(x, digits=digits, scientific=scientific, ...), ...)

  invisible(x)
}
