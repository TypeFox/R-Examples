.df <- function(x) attr(logLik(x), "df")

#' Returns degrees of freedom, log likelihood, R-squared, AIC, BIC and
#' adjusted R-squared.
#'
#' @param object ensemble of models
#' @param ... other arguments ignored
#' @keywords regression
#' @export
summary.ensemble <- function(object, ...) {
  fits <- data.frame(t(sapply(object, function(mod) {
    sum <- summary(mod)
    c(
      df = .df(mod),
      logL = logLik(mod),
      AIC = -AIC(mod),
      BIC = -AIC(mod, k=log(length(fitted(mod)))),
      R2 = sum$r.squared,
      adjR2 = sum$adj.r.squared,
      n = length(mod$residuals)
    )
  })))
  fits$model <- factor(names(object))
  rownames(fits) <- paste("m", fits$model, sep="")
  fits
}

#' Calculcate coefficients for all models in ensemble.
#' Returns raw, t-value, absolute t-value, and standardised coefficent values.
#'
#' @param object ensemble of models
#' @param ... other arguments ignored
#' @keywords regression
#' @export
coef.ensemble <- function(object, ...) {
  coefs <- ldply(object, coef_simple, data = attr(object, "data"))
  names(coefs)[1] <- "model"
  coefs$model <- factor(coefs$model)

  all <- expand.grid(
      model = unique(coefs$model),
      variable = unique(coefs$variable))
  coefs <- join(all, coefs, by = c("model", "variable"))
  coefs[is.na(coefs)] <- 0
  rownames(coefs) <- paste("m", coefs$model, "v", as.numeric(coefs$variable),
    sep = "")
  class(coefs) <- c("variable_ensemble", class(coefs))

  coefs
}

coef_simple <- function(model, data) {
  trunc <- function(x, trunc) ifelse(abs(x) > trunc, sign(x) * trunc, x)

  coefs <- data.frame(
    names(coef(model))[-1],
    summary(model)$coefficients[-1, c(1, 3), drop=FALSE]
  )
  names(coefs) <- c("variable", "raw", "t")
  transform(coefs,
    t = trunc(t,3),
    abst = abs(trunc(t, 3)),
    std = stdcoef(model, data)[-1]
  )
}

#' Summarise variable ensemble.
#'
#' Provides variable level statistics.
#'
#' @param object ensemble of models
#' @param ... other arguments ignored
#' @keywords regression
#' @export
summary.variable_ensemble <- function(object, ...) {
  coefs <- subset(object, raw != 0)

  ddply(coefs, "variable", summarise,
     raw_mean = mean(raw),
     raw_sd = sd(raw),
     t_mean = mean(raw),
     t_sd = sd(t),
     std_mean = mean(std),
     std_sd = sd(std),
     n = length(model))
}
globalVariables(c("std", "model"))

# Calculcate standardised coefficients for a model
stdcoef <- function(model, data = model$model) {
  data[] <- lapply(data, scale)
  coef(update(model, . ~ ., data = data))
}

#' Calculate residuals for all models in ensemble.
#'
#' @return data.frame of class \code{resid_ensemble}
#' @seealso \code{\link{summary.resid_ensemble}}
#' @param object ensemble of models
#' @param ... other arguments ignored
#' @keywords regression
#' @export
residuals.ensemble <- function(object, ...) {
  resids <- do.call(rbind, mapply(function(mod, name) {
    data.frame(
      obs=gsub("\\.[0-9]+$", "", names(resid(mod))),
      resid=resid(mod),
      rstudent = rstudent(mod),
      fitted=fitted(mod),
      true=fitted(mod) + resid(mod),
      influence.measures(mod)$infmat[, c("dffit", "cov.r", "cook.d", "hat")],
      model=name)
  }, object, names(object), SIMPLIFY=FALSE))
  resids$model <- factor(resids$model)
  rownames(resids) <- 1:nrow(resids)

  scores <- tapply(resids$resid, resids$obs, mean)
  resids$obs <- factor(resids$obs, levels = names(scores)[order(scores)])
  class(resids) <- c("resid_ensemble", class(resids))
  attr(resids, "data") <- attr(object, "data")
  resids
}

#' Summarise residuals from ensemble.
#'
#' @param object model residuals from \code{\link{residuals.ensemble}}
#' @param data associated data set
#' @param ... other arguments ignored
#' @keywords regression
#' @export
summary.resid_ensemble <- function(object, data = attr(object, "data"), ...) {
  s <- ddply(object, "obs", summarise,
    mean = mean(rstudent),
    sd = sd(rstudent),
    n = length(rstudent))

  if (!is.null(data)) {
    data$obs <- rownames(data)
    rownames(data) <- NULL
    s <- join(s, data, by = "obs")
  }
  s
}
