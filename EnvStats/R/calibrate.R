calibrate <-
function (formula, data, test.higher.orders = TRUE, max.order = 4, 
    p.crit = 0.05, F.test = "partial", weights, subset, na.action, 
    method = "qr", model = FALSE, x = FALSE, y = FALSE, contrasts = NULL, 
    warn = TRUE, ...) 
{
    fun.call <- match.call()
    fun.call[[1]] <- as.name("lm")
    fun.args <- is.element(arg.names(fun.call), names(formals(lm)))
    fun.call <- fun.call[c(TRUE, fun.args)]
    fit <- eval(fun.call)
    name.pred <- attr(terms(fit), "term.labels")
    E <- length(name.pred)
    if (!test.higher.orders) {
        if (E >= 2) {
            if (length(grep(name.pred[1], name.pred)) < E) {
                stop(paste("All predictor variables in the model must be functions of", 
                  "a single variable; for example, x, x^2, etc."))
            }
            pred <- model.frame(fit)[, name.pred[1]]
            dc.pred <- data.class(pred)
            if (!((dc.pred == "AsIs" & is.numeric(pred)) || dc.pred == 
                "numeric")) {
                stop("The single variable that all predictors are functions of must be numeric.")
            }
        }
        else {
            pred <- model.frame(fit)[, name.pred]
            dc.pred <- data.class(pred)
            if (!((dc.pred == "AsIs" & is.numeric(pred)) || dc.pred == 
                "numeric")) {
                stop("The single predictor variable must be numeric.")
            }
        }
        if (!x) {
            fit <- update(fit, x = TRUE)
        }
        return(fit)
    }
    if (E != 1) {
        stop(paste("When test.higher.orders=TRUE, there can be only one", 
            "predictor variable in the initial calibration model"))
    }
    pred <- model.frame(fit)[, name.pred]
    dc.pred <- data.class(pred)
    if (!((dc.pred == "AsIs" & is.numeric(pred)) || dc.pred == 
        "numeric")) {
        stop("The single predictor variable must be numeric.")
    }
    if (!is.numeric(max.order) || length(max.order) != 1 || max.order < 
        1 || max.order != trunc(max.order)) {
        stop("The argument 'max.order' must be an integer greater than 0")
    }
    if (!is.numeric(p.crit) || length(p.crit) != 1 || p.crit <= 
        0 || p.crit >= 1) {
        stop("The argument 'p.crit' must be a numeric scalar greater than 0 and less than 1")
    }
    F.test <- match.arg(F.test, c("lof", "partial"))
    n.pred <- length(pred)
    n.pred.unique <- length(unique(pred))
    max.order.new <- min(max.order, n.pred.unique - 1)
    if (max.order.new < max.order) {
        if (warn) {
            warning(paste("The argument 'max.order' was reset from", 
                max.order, "to", max.order.new, "because there are only", 
                n.pred.unique, "unique values of the single predictor variable."))
        }
        max.order <- max.order.new
    }
    if (E < max.order) {
        if (F.test == "lof" && (n.pred.unique < n.pred) && fit$df.residual > 
            (n.pred - n.pred.unique) + 1) {
            aov.table <- anovaPE(fit)
            index <- grep("Lack of Fit", row.names(aov.table))
            lof.p <- aov.table[index, "Pr(>F)"]
            try.new <- lof.p < p.crit
            while (try.new) {
                E <- E + 1
                formula.new <- eval(parse(text = paste(". ~  . + I(", 
                  name.pred, "^", eval(E), ")", sep = "", collapse = "")))
                fit.new <- update(fit, formula = formula.new)
                if (any(is.na(fit.new$coef))) {
                  if (warn) {
                    warning(paste("Final model of order", E - 
                      1, "because of singularities in higher order models."))
                  }
                  break
                }
                if (fit.new$df.residual > (n.pred - n.pred.unique)) {
                  fit <- fit.new
                  aov.table <- anovaPE(fit)
                  index <- grep("Lack of Fit", row.names(aov.table))
                  lof.p <- aov.table[index, "Pr(>F)"]
                  try.new <- (E < max.order) && (lof.p < p.crit)
                }
                else {
                  partial.F.p <- anova(fit, fit.new)[2, "Pr(>F)"]
                  L1 <- partial.F.p < p.crit
                  if (L1) {
                    fit <- fit.new
                  }
                  try.new <- L1 && (E < max.order)
                }
            }
        }
        else {
            try.new <- TRUE
            while (try.new) {
                E <- E + 1
                formula.new <- eval(parse(text = paste(". ~  . + I(", 
                  name.pred, "^", eval(E), ")", sep = "", collapse = "")))
                fit.new <- update(fit, formula = formula.new)
                if (any(is.na(fit.new$coef))) {
                  if (warn) {
                    warning(paste("Final model of order", E - 
                      1, "because of singularities in higher order models."))
                  }
                  break
                }
                partial.F.p <- anova(fit, fit.new)[2, "Pr(>F)"]
                L1 <- partial.F.p < p.crit
                if (L1) {
                  fit <- fit.new
                }
                try.new <- L1 && (E < max.order)
            }
        }
    }
    if (!x) 
        fit <- update(fit, x = TRUE)
    class(fit) <- c("calibrate", "lm")
    fit
}
