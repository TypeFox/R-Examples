intcox <-
function (formula = formula(data), data = parent.frame(), subset,
    na.action, x = FALSE, y = TRUE, epsilon = 1e-04, itermax = 10000,
    no.warnings = FALSE)
{
#
#   preparation
#
    copy.data <- data
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "copy.data", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    m[[1]] <- as.name("model.frame")
    Terms <- if (missing(copy.data))
        terms(formula)
    else terms(formula, data = copy.data)
    m$formula <- Terms
    m <- eval(m, parent.frame())
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    attr(Terms, "intercept") <- 1
    type <- attr(Y, "type")
    if (type != "interval")
        stop("Invalid survival type (only interval censored data)")
    else {
        if (any(Y[, 3] != 3 && Y[, 3] != 0))
            stop("Invalid cens status")
        else {
            if (any(Y[, 3] == 3))
                Y <- cbind(Y[, 1:2], Y[, 3])
            else Y <- cbind(Y[, 1], Y[, 3])
        }
    }
    copy.data$left <- Y[, 1]
    copy.data$right <- Y[, 2]
    copy.data$cens <- Y[, 3]
    copy.data$mix <- ifelse(copy.data$cens == 3, copy.data$right,
        copy.data$left)                             # start with the right ends if observed
    ord <- order(copy.data$mix, 3 - copy.data$cens) #sort, first by right interval limit, then by censoring, that the Breslow-estimator can work
    copy.data <- copy.data[ord, ]
    X <- model.matrix(Terms, m)
    X <- X[, -1, drop = FALSE]
    covar <- X[ord, , drop = FALSE]
#
#   calculation
#
    if (length(attr(Terms, "term.label")) == 0)
        stop("Invalid type of model")
    else intcox.new <- intcox.fit(formula, copy.data, covar, epsilon,
        itermax)
#
#   output
#
    if (is.character(intcox.new)) {
        intcox.new <- list(fail = intcox.new)
        class(intcox.new) <- "coxph"
    }
    else {
        intcox.new$n <- nrow(Y)
        class(intcox.new) <- "coxph"
        intcox.new$terms <- Terms
        intcox.new$assign <- assign
        intcox.new$var <- matrix(NA, ncol = 1)
        intcox.new$wald.test <- NA
        intcox.new$score <- NA
        intcox.new$loglik <- intcox.new$likeli.vec[length(intcox.new$likeli.vec)]
        intcox.new$linear.predictors <- as.vector(intcox.new$coefficients *
            t(covar))
        intcox.new$residuals <- NA
        intcox.new$means <- apply(covar, 2, mean)
        na.action <- attr(m, "na.action")
        if (length(na.action))
            intcox.new$na.action <- na.action
        if (x) {
            intcox.new$x <- covar
        }
        if (y)
            intcox.new$y <- Y
    }
    intcox.new$formula <- formula(Terms)
    intcox.new$call <- call
    intcox.new$method <- NA
    if (no.warnings==FALSE) {
        if (intcox.new$termination == 4)
            cat("inside precondition(s) at iteration = ", intcox.new$iter,
                "not fulfilled\n")
        if (intcox.new$termination == 2)
            cat("no improvement of likelihood possible, iteration = ",
                intcox.new$iter, "\n")
        if (intcox.new$termination == 3)
            cat("algorithm did not converge - maximum number of iteration reached, itermax = ",
                intcox.new$iter, "\n")
    }
    return(intcox.new)
}
