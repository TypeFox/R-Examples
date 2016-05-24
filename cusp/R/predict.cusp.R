as.lm.cusp <-
function (object, what = c("y", "alpha", "beta"))
  {
      what = match.arg(what)
      object$coefficients = object$coefficients[grep(sub("a", switch(what,
          alpha = "a", beta = "b", y = "w"), "^a\\["), names(object$coefficients))]
      object$terms = object$terms[[what]]
      object$call$formula = object$call[[switch(what, alpha = "alpha",
          beta = "beta", y = "formula")]]
      object$model = object$model[[what]]
      object$rank = attr(object$rank, "ranks")[switch(what, alpha = "ranka",
          beta = "rankb", y = "ranky")]
      object$contrasts = object$contrasts[[what]]
      object$xlevels = object$xlevels[[what]]
      class(object) = "lm"
      object
  }


predict.cusp <-
function (object, newdata, se.fit = FALSE, interval = c("none",
    "confidence", "prediction"), level = 0.95, type = c("response",
    "terms"), terms = NULL, na.action = na.pass, pred.var = res.var/weights,
    weights = 1, method = c("delay", "maxwell", "expected"),
    keep.linear.predictors = FALSE, ...)
    {
        if (!all(missing(se.fit), missing(interval), missing(level),
            missing(type), missing(terms), missing(pred.var), missing(weights)))
            .NotYetUsed(intersect(names(match.call()), c("se.fit",
                "interval", "level", "type", "terms", "pred.var",
                "weights")), error = TRUE)
        if (se.fit || interval != "none") {
          res.var <- (1-summary(object)$r2cusp.r.squared[1]) * var(drop(object$y))
        }
        method = match.arg(method)
        if (missing(newdata)) {
            alpha = object$linear.predictors[, "alpha"]
            beta = object$linear.predictors[, "beta"]
        }
        else {
            alpha = predict(as.lm.cusp(object, "alpha"), newdata,
                se.fit = se.fit, interval = interval, level = level,
                type = type, terms = terms, na.action = na.action,
                pred.var = pred.var, weights = weights, ...)
            alpha = if (is.list(alpha))
                alpha$fit
            else as.matrix(alpha)[, 1]
            beta = predict(as.lm.cusp(object, "beta"), newdata, se.fit = se.fit,
                interval = interval, level = level, type = type,
                terms = terms, na.action = na.action, pred.var = pred.var,
                weights = weights, ...)
            beta = if (is.list(beta))
                beta$fit
            else as.matrix(beta)[, 1]
        }
        if (method == "delay") {
            if (missing(newdata)) {
                y = drop(object$y)
            }
            else {
                y = predict(as.lm.cusp(object, "y"), newdata, se.fit = se.fit,
                    interval = interval, level = level, type = type,
                    terms = terms, na.action = na.action, pred.var = pred.var,
                    weights = weights, ...)
                y = if (is.list(y))
                    y$fit
                else as.matrix(y)[, 1]
            }
            pred = cusp.expected(alpha, beta, y, method = method)
        }
        else {
            pred = cusp.expected(alpha, beta, method = method)
        }
        if (keep.linear.predictors && !missing(newdata)) {
            if (method == "delay")
                attr(pred, "data") = cbind(newdata, alpha = alpha,
                    beta = beta, y = y)
            else attr(pred, "data") = cbind(newdata, alpha = alpha,
                beta = beta)
        }
        pred
    }
