.check.alpha <-
function (call, alpha, n.stud) 
{
    if (!is.numeric(alpha)) {
        .stop(call, "alpha must be a numeric vector")
    }
    if (any(is.na(alpha))) {
        .stop(call, "alpha cannot have missing values")
    }
    if (any(alpha <= 0) || any(alpha >= 1)) {
        .stop(call, "alpha cannot be <= 0 or >= 1")
    }
    if (length(alpha) == 1) {
        return(rep(alpha, n.stud))
    }
    if (length(alpha) != n.stud) {
        .stop(call, "alpha has an incorrect length")
    }
    alpha
}
.check.data <-
function (call, data, n.stud) 
{
    if (is.vector(data)) {
        if (length(data) != n.stud) {
            .stop(call, "data has an incorrect length")
        }
        if (any(is.na(data))) {
            .stop(call, "data cannot have missing values. Impute missing values (using for example the R package 'mi'), call 'meta' for each imputation, and combine all imputations")
        }
        return(model.matrix(as.formula("~data")))
    }
    else if (!is.matrix(data) && !is.data.frame(data)) {
        .stop(call, "data must be a vector, matrix or data.frame")
    }
    if (!nrow(data)) {
        return(matrix(1, n.stud))
    }
    if (nrow(data) != n.stud) {
        .stop(call, "data has an incorrect number of rows")
    }
    if (any(is.na(data))) {
        .stop(call, "data cannot have missing values. Impute missing values (using for example the R package 'mi'), call 'meta' for each imputation, and combine all imputations")
    }
    xnames <- paste0("data", 1:ncol(data))
    colnames(data) <- xnames
    model.matrix(as.formula(paste0("~", paste(xnames, collapse = "+"))), 
        data)
}
.check.formula <-
function (call, formula, n.stud) 
{
    if (!inherits(formula, "formula")) {
        .stop(call, "formula must be a formula")
    }
    terms <- terms(formula)
    xnames <- attr(terms, "term.labels")
    if (length(xnames) == 0) {
        return(list(formula = "~ 1", n.coef = 1, matrix = matrix(1, 
            n.stud), labels = "(Mean)"))
    }
    formula <- paste("~", paste(xnames, collapse = " + "))
    if (!attr(terms, "intercept")) {
        warning("You have specified a regression though the origin")
        formula <- paste(formula, "- 1")
    }
    X <- model.matrix(as.formula(formula), parent.frame(2))
    if (nrow(X) != n.stud) {
        .stop(call, "Independent variables of the formula have an incorrect length")
    }
    if (any(is.na(X))) {
        .stop(call, "Independent variables of the formula cannot have missing values. Impute missing values (using for example the R package 'mi'), call 'meta' for each imputation, and combine all imputations")
    }
    list(formula = formula, n.coef = ncol(X), matrix = matrix(c(X), 
        n.stud), labels = colnames(X))
}
.check.hypotheses <-
function (call, hypotheses, model) 
{
    n.coef <- model$n.coef
    if (is.null(hypotheses)) {
        matrixs <- list()
        for (i in 1:n.coef) {
            matrixs[[i]] <- matrix(c(rep(0, i - 1), 1, rep(0, 
                n.coef - i)), 1)
        }
        if (n.coef > 2) {
            n.hyp <- n.coef + 1
            labels <- c(model$labels, "omnibus")
            matrixs[[n.coef + 1]] <- matrix(c(0, rep(1, n.coef - 
                1)), 1)
        }
        else {
            n.hyp <- n.coef
            labels <- model$labels
        }
        return(list(n.hyp = n.hyp, labels = labels, matrixs = matrixs))
    }
    else if (is.list(hypotheses)) {
        n.hyp <- length(hypotheses)
        if (!n.hyp) {
            .stop(call, "No hypotheses defined")
        }
        labels <- names(hypotheses)
        matrixs <- list()
        for (i in 1:n.hyp) {
            hypothesis <- hypotheses[[i]]
            if (is.matrix(hypothesis)) {
                if (ncol(hypothesis) != n.coef) {
                  .stop(call, "Wrong number of columns")
                }
            }
            else if (is.numeric(hypothesis)) {
                if (length(hypothesis) != n.coef) {
                  .stop(call, "Wrong vector length")
                }
                hypothesis <- matrix(hypothesis, 1)
            }
            else {
                .stop(call, "Numeric vector or matrix expected")
            }
            matrixs[[i]] <- hypothesis
        }
    }
    else {
        .stop(call, "hypotheses should be a list")
    }
    list(n.hyp = n.hyp, labels = labels, matrixs = matrixs)
}
.check.labels <-
function (call, labels, n.stud) 
{
    if (is.factor(labels)) {
        labels <- as.character(labels)
    }
    else if (!is.vector(labels)) {
        .stop(call, "labels must be a vector")
    }
    if (length(labels) == 1) {
        return(paste0(labels, 1:n.stud))
    }
    if (length(labels) != n.stud) {
        .stop(call, "labels has an incorrect length")
    }
    labels
}
.check.n <-
function (call, n, min.n, n.stud) 
{
    if (!is.numeric(n)) {
        .stop(call, "n must be a numeric vector")
    }
    if (any(is.na(n))) {
        .stop(call, "n cannot have missing values")
    }
    if (any(n < min.n)) {
        .stop(call, paste("n cannot be <", min.n))
    }
    if (length(n) == 1) {
        return(rep(n, n.stud))
    }
    if (length(n) != n.stud) {
        .stop(call, "n has an incorrect length")
    }
    n
}
coef.meta.nsue <-
function (object, ...) 
{
    coef <- object$model$coef
    se <- object$model$se
    table <- cbind(coef, se, object$model$z, object$model$p.value, 
        coef + qnorm(0.025) * se, coef + qnorm(0.975) * se)
    colnames(table) <- c("Estimate", "Std. Error", "z value", 
        "Pr(>|z|)", "CI(low)", "CI(up)")
    table
}
.d_j <-
function (x) 
{
    j <- gamma(x/2)/(sqrt(x/2) * gamma((x - 1)/2))
    na.j <- which(is.na(j))
    j[na.j] <- 1 - 3/(4 * x[na.j] - 1)
    j
}
.elliptic.q <-
function (x, y, p = 0.95, col = "#cccccc", segments = 51) 
{
    center <- c(mean(x), mean(y))
    shape <- cov(cbind(x, y))
    radius <- sqrt(2 * qf(p, 2, length(x) - 1))
    angles <- (0:segments) * 2 * pi/segments
    circle <- cbind(cos(angles), sin(angles))
    choleski <- chol(shape, pivot = TRUE)
    polygon(t(center + radius * t(circle %*% choleski[, order(attr(choleski, 
        "pivot"))])), col = col, border = NA)
}
fitted.meta.nsue <-
function (object, ...) 
{
    object$model$matrix %*% object$model$coef
}
forest <-
function (x, ...) 
UseMethod("forest")
forest.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    if (x$rm$n.stud < x$known$n.stud + x$unknown$n.stud) {
        .warning("This plot shows repeated measures as separate studies")
    }
    measure = x$measure
    known = x$known$i
    unknown = x$unknown$i
    unknown.n.stud <- x$unknown$n.stud
    n.hyp <- x$hypotheses$n.hyp
    n.stud <- x$known$n.stud + unknown.n.stud
    n <- n.stud + n.hyp
    labels <- c(x$labels[unknown], x$labels[known], x$hypotheses$labels)
    pos.y <- c(n.stud + 2 - c(unknown, known), 0:(1 - n.hyp))
    if (x$unknown$n.stud) {
        y <- apply(x$unknown$y, 1, mean)
        y.low <- apply(x$unknown$y, 1, function(x) {
            quantile(x, 0.025)
        })
        y.upp <- apply(x$unknown$y, 1, function(x) {
            quantile(x, 0.975)
        })
        if (measure == "cor" || measure == "cor in smd") {
            y.se <- sqrt(x$y.var[unknown])
        }
        if (measure == "smc" || measure == "smd") {
            y.se <- sqrt(.pool.var(x$unknown$y, x$y2var_k1[unknown] + 
                x$y2var_k2[unknown] * x$unknown$y^2))
        }
    }
    else {
        y <- y.low <- y.upp <- y.se <- c()
    }
    y <- c(y, x$known$y, x$hypotheses$coef)
    if (measure == "cor" || measure == "cor in smd") {
        y.se <- c(y.se, sqrt(x$y.var[known]), x$hypotheses$se)
    }
    if (measure == "smc" || measure == "smd") {
        y.se <- c(y.se, sqrt(x$y2var_k1[known] + x$y2var_k2[known] * 
            x$known$y^2), x$hypotheses$se)
    }
    ci.low <- y + qnorm(0.025) * y.se
    ci.upp <- y + qnorm(0.975) * y.se
    if (measure == "cor" || measure == "cor in smd") {
        y <- tanh(y)
        y.low <- tanh(y.low)
        y.upp <- tanh(y.upp)
        ci.low <- tanh(ci.low)
        ci.upp <- tanh(ci.upp)
    }
    lwd <- 1/y.se
    lwd <- sqrt(9 + 216 * (lwd - min(lwd))/(max(lwd) - min(lwd)))
    ci.text <- paste0(.format.2(y), " [ ", .format.2(ci.low), 
        ", ", .format.2(ci.upp), " ] ", .format.sign(2 * pnorm(-abs(y/y.se))))
    plot.new()
    xlim <- c(-2.5 - max(strwidth(labels, units = "inches")), 
        max(strwidth(ci.text, units = "inches")) + 2.5)
    ylim <- c(-1 - n.hyp, n)
    plot.window(xlim = xlim, ylim = ylim)
    xthr <- .signif.up(max(abs(c(quantile(ci.low, 0.1), quantile(ci.upp, 
        0.9)))), 1)
    lines(rep(0, 2), c(n + 0.5, -n.hyp), col = "#bbbbbb", lty = 1)
    lines(c(-2, 2), rep(-n.hyp, 2), col = "#bbbbbb", lty = 1)
    for (pos.x in -2:2) {
        lines(rep(pos.x, 2), c(-n.hyp, -n.hyp - 0.3), col = "#bbbbbb", 
            lty = 1)
        text(pos.x, -n.hyp - 0.5, .format.2(pos.x/2 * xthr), 
            pos = 1, col = "#bbbbbb")
    }
    for (i in 1:n) {
        pos.y_i <- pos.y[i]
        y_i <- y[i]
        ci.low_i <- ci.low[i]
        ci.upp_i <- ci.upp[i]
        if (i > unknown.n.stud) {
            col <- "#000000"
        }
        else {
            y.low_i <- y.low[i]
            y.upp_i <- y.upp[i]
            if (y.upp_i > -xthr && y.low_i < xthr) {
                lines(c(max(y.low_i/xthr * 2, -2), min(y.upp_i/xthr * 
                  2, 2)), rep(pos.y_i, 2), lwd = lwd[i], col = "#dddddd")
            }
            col <- "#a7a7a7"
        }
        if (y_i > -xthr && y_i < xthr) {
            lines(rep(y_i/xthr * 2, 2), rep(pos.y_i, 2), lwd = lwd[i], 
                col = col)
        }
        if (ci.upp_i > -xthr && ci.low_i < xthr) {
            lines(c(max(ci.low_i/xthr * 2, -2), min(ci.upp_i/xthr * 
                2, 2)), rep(pos.y_i, 2), lend = 2, col = col)
            if (ci.low_i > -xthr) {
                lines(rep(ci.low_i/xthr * 2, 2), pos.y_i + c(-0.15, 
                  0.15), lend = 2, col = col)
            }
            if (ci.upp_i < xthr) {
                lines(rep(ci.upp_i/xthr * 2, 2), pos.y_i + c(-0.15, 
                  0.15), lend = 2, col = col)
            }
        }
        text(-2.1, pos.y_i, labels[i], pos = 2, col = col)
        text(2.1, pos.y_i, ci.text[i], pos = 4, col = col)
    }
    width <- round(diff(xlim))
    height <- round(diff(ylim)/3)
    cat("\n")
    cat("Use pdf(filename, width, height) before calling forest to save it.\n")
    cat("The optimal width and height of this plot is ~", width, 
        " x ~", height, " inches.\n", sep = "")
    cat("\n")
    invisible(list(optimal.width = width, optimal.height = height))
}
.format.0pos <-
function (x) 
{
    formatC(x, 0, width = 2, format = "f")
}
.format.1 <-
function (x) 
{
    formatC(x, 1, width = 4, format = "f")
}
.format.1pos <-
function (x) 
{
    formatC(x, 1, width = 3, format = "f")
}
.format.2 <-
function (x) 
{
    formatC(x, 2, width = 5, format = "f")
}
.format.2pos <-
function (x) 
{
    formatC(x, 2, width = 4, format = "f")
}
.format.3 <-
function (x) 
{
    formatC(x, 3, width = 6, format = "f")
}
.format.4 <-
function (x) 
{
    formatC(x, 4, width = 7, format = "f")
}
.format.4pos <-
function (x) 
{
    formatC(x, 4, width = 6, format = "f")
}
.format.measure <-
function (measure) 
{
    switch(measure, cor = "correlation", smc = "standardized mean change", 
        smd = "standardized mean difference")
}
.format.perc2 <-
function (x) 
{
    formatC(100 * x, 2, width = 5, format = "f")
}
.format.prob <-
function (p) 
{
    p <- formatC(p, digits = 4, format = "f")
    p[p == "0.0000"] <- "<.0001"
    p
}
.format.sign <-
function (p) 
{
    symnum(p, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
        "**", "*", ".", " "), na = FALSE, corr = FALSE)
}
funnel <-
function (x, ...) 
UseMethod("funnel")
funnel.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    if (x$rm$n.stud < x$known$n.stud + x$unknown$n.stud) {
        .warning("This analysis does not take repeated measures into account")
    }
    measure = x$measure
    known <- x$known$i
    unknown <- x$unknown$i
    fitted <- fitted(x)
    known.res <- x$known$y - fitted[known]
    if (measure == "cor" || measure == "cor in smd") {
        known.se <- sqrt(x$y.var[known])
    }
    if (measure == "smc" || measure == "smd") {
        known.se <- sqrt(x$y2var_k1[known] + x$y2var_k2[known] * 
            x$known$y^2)
    }
    if (x$unknown$n.stud) {
        unknown.fitted <- fitted[unknown]
        unknown.res <- apply(x$unknown$y, 1, mean) - unknown.fitted
        if (measure == "cor" || measure == "cor in smd") {
            unknown.se <- sqrt(x$y.var[unknown])
        }
        if (measure == "smc" || measure == "smd") {
            unknown.se <- sqrt(apply(x$y2var_k1[unknown] + x$y2var_k2[unknown] * 
                x$unknown$y^2, 1, mean))
        }
        max.se <- .signif.up(max(c(known.se, unknown.se)), 1)
    }
    else {
        max.se <- .signif.up(max(known.se), 1)
    }
    ci <- qnorm(0.975) * max.se
    plot(NA, NA, type = "n", xlim = 1.3 * c(-ci, ci), ylim = c(max.se, 
        0), lty = 2, frame.plot = FALSE, xlab = "Residual effect size", 
        ylab = "Standard error")
    ci.x <- c(-ci, 0, ci)
    ci.y <- c(max.se, 0, max.se)
    polygon(c(ci.x, rep(1.3 * ci, 2), rep(-1.3 * ci, 2)), c(ci.y, 
        max.se, 0, 0, max.se), col = "#fcfcfc", border = "#dddddd")
    if (x$unknown$n.stud) {
        for (i in 1:x$unknown$n.stud) {
            if (measure == "cor" || measure == "cor in smd") {
                lines(c(max(quantile(x$unknown$y[i, ] - unknown.fitted[i], 
                  0.025), -1.3 * ci), min(quantile(x$unknown$y[i, 
                  ] - unknown.fitted[i], 0.975), 1.3 * ci)), 
                  rep(sqrt(x$y.var[unknown][i]), 2), lwd = 7, 
                  col = "#dddddd")
            }
            if (measure == "smc" || measure == "smd") {
                .elliptic.q(x$unknown$y[i, ] - unknown.fitted[i], 
                  sqrt(x$y2var_k1[unknown][i] + x$y2var_k2[unknown][i] * 
                    x$unknown$y[i, ]^2), col = "#dddddd")
            }
        }
    }
    lines(ci.x, ci.y, lty = 2)
    lines(c(0, 0), c(max.se, 0), lty = 2)
    if (x$unknown$n.stud) {
        for (i in 1:x$unknown$n.stud) {
            lines(rep(unknown.res[i], 2), rep(unknown.se[i], 
                2), lwd = 7, col = "#a7a7a7")
        }
    }
    for (i in 1:x$known$n.stud) {
        lines(rep(known.res[i], 2), rep(known.se[i], 2), lwd = 7, 
            col = "#000000")
    }
    cat("\n")
    cat("Use pdf(filename) before calling funnel to save it.\n")
    cat("\n")
}
leave1out <-
function (x, ...) 
UseMethod("leave1out")
leave1out.nsue <-
function (x, data = data.frame(), formula = ~1, hypotheses = NULL, 
    n.imp = 50, n.bins = 200, maxiter = 200, tol = 1e-06, ...) 
{
    call <- match.call()
    y <- x$y
    measure <- x$measure
    n.stud <- length(y)
    X <- .check.data(call, data, n.stud)
    model <- .check.formula(call, formula, n.stud)
    hypotheses <- .check.hypotheses(call, hypotheses, model)
    if (n.imp < 2) {
        .stop(call, "At least 2 imputations XXX")
    }
    if (length(unique(x$labels)) < n.stud) {
        .warning("This analysis understand repeated measures as separate studies")
    }
    nsue_i <- x
    model_i <- model
    obj <- list()
    for (i in 1:n.stud) {
        nsue_i$y <- x$y[-i]
        nsue_i$y_lo <- x$y_lo[-i]
        nsue_i$y_up <- x$y_up[-i]
        if (measure == "cor" || measure == "cor in smd") {
            nsue_i$y.var <- x$y.var[-i]
        }
        if (measure == "smc" || measure == "smd") {
            nsue_i$y2var_k1 <- x$y2var_k1[-i]
            nsue_i$y2var_k2 <- x$y2var_k2[-i]
        }
        if (measure == "cor in smd") {
            nsue_i$smd = x$smd[-i, ]
        }
        nsue_i$labels <- x$labels[-i]
        class(nsue_i) <- "nsue"
        model_i$matrix <- as.matrix(model$matrix[-i, ])
        obj[[i]] <- list(study = x$labels[i], meta.nsue = .meta.nsue(nsue_i, 
            as.matrix(X[-i, ]), model_i, hypotheses, n.imp, n.bins, 
            maxiter, tol))
    }
    class(obj) <- "leave1out.nsue"
    obj
}
linearHypothesis <-
function (x, ...) 
UseMethod("linearHypothesis")
.linearHypothesis.meta.nsue <-
function (x) 
{
    n.hyp <- x$hypotheses$n.hyp
    coef <- c()
    se <- c()
    z <- c()
    chisq <- c()
    df <- c()
    p.value <- c()
    for (i in 1:n.hyp) {
        matrix <- x$hypotheses$matrixs[[i]]
        coef_i <- c(matrix %*% x$model$coef)
        cov <- matrix %*% x$model$cov %*% t(matrix)
        if (nrow(matrix) == 1) {
            se_i <- sqrt(cov)
            z_i <- coef_i/se_i
            coef[i] <- coef_i
            se[i] <- se_i
            z[i] <- z_i
            chisq[i] <- NA
            df[i] <- NA
            p.value[i] <- 2 * pnorm(-abs(z_i))
        }
        else {
            qr <- qr(cov)
            pivot <- qr$pivot[1:qr$rank]
            chisq_i <- c(t(coef_i[pivot]) %*% solve(cov[pivot, 
                pivot]) %*% coef_i[pivot])
            df_i <- length(pivot)
            coef[i] <- NA
            se[i] <- NA
            z[i] <- NA
            chisq[i] <- chisq_i
            df[i] <- df_i
            p.value[i] <- 1 - pchisq(chisq_i, df_i)
        }
    }
    x$hypotheses$coef <- coef
    x$hypotheses$se <- se
    x$hypotheses$z <- z
    x$hypotheses$chisq <- chisq
    x$hypotheses$df <- df
    x$hypotheses$p.value <- p.value
    x
}
linearHypothesis.meta.nsue <-
function (x, hypotheses, ...) 
{
    call <- match.call()
    if (!inherits(x, "meta.nsue")) {
        .stop(call, "The argument must be a 'meta.nsue' object")
    }
    x$hypotheses <- .check.hypotheses(call, hypotheses, x$model)
    .linearHypothesis.meta.nsue(x)
}
meta <-
function (x, ...) 
UseMethod("meta")
metabias <-
function (x, ...) 
UseMethod("metabias")
metabias.meta.nsue <-
function (x, maxiter = 100, tol = 1e-06, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    if (x$rm$n.stud < x$known$n.stud + x$unknown$n.stud) {
        .warning("This analysis does not take repeated measures into account")
    }
    measure <- x$measure
    known <- x$known$i
    unknown <- x$unknown$i
    X <- x$model$matrix
    n.coef_j <- x$model$n.coef + 1
    se_coef <- c()
    se_coef_var <- c()
    for (j in 1:x$unknown$n.imp) {
        y <- c(x$known$y, x$unknown$y[, j])
        if (measure == "cor" || measure == "cor in smd") {
            y.var <- c(x$known$y.var, x$unknown$y.var[, j])
        }
        if (measure == "smc" || measure == "smd") {
            y.var <- c(x$y2var_k1[known] + x$y2var_k2[known] * 
                x$known$y^2, x$y2var_k1[unknown] + x$y2var_k2[unknown] * 
                x$unknown$y[, j]^2)
        }
        X_j <- cbind(X, sqrt(y.var))
        W <- diag(1/(y.var + .tau2.reml(x$heterogeneity$tau2, 
            y, y.var, X_j, maxiter, tol)))
        inv_XtWX <- solve(t(X_j) %*% W %*% X_j)
        se_coef <- c(se_coef, (inv_XtWX %*% t(X_j) %*% W %*% 
            y)[n.coef_j])
        se_coef_var <- c(se_coef_var, diag(inv_XtWX)[n.coef_j])
    }
    z <- mean(se_coef)/sqrt(.pool.var(se_coef, se_coef_var))
    names(z) <- "z"
    p <- 2 * pnorm(-abs(z))
    x <- list(method = "'meta.nsue' regression test for funnel plot asymmetry", 
        data.name = as.character(match.call()[2]), statistic = z, 
        p.value = p)
    class(x) <- "htest"
    x
}
metalm <-
function (x, ...) 
UseMethod("metalm")
.metalm.meta.nsue <-
function (x, maxiter, tol) 
{
    measure <- x$measure
    known <- x$known$i
    unknown <- x$unknown$i
    n.imp <- x$unknown$n.imp
    rm <- x$rm
    X <- x$model$matrix
    n.coef <- x$model$n.coef
    df <- rm$n.stud - n.coef
    tau2 <- x$heterogeneity$tau2
    rm.M <- t(apply(rm$M, 1, function(x) {
        x/sum(x)
    }))
    aX <- rm.M %*% X
    y <- c()
    y[known] <- x$known$y
    if (measure == "cor" || measure == "cor in smd") {
        y.var <- x$y.var/rm$weights
        ay.var <- apply(rm$M, 1, function(xx) {
            mean(y.var[which(xx == 1)])/sum(xx)
        })
    }
    if (measure == "smc" || measure == "smd") {
        y2var_k1 <- x$y2var_k1/rm$weights
        ay2var_k1 <- apply(rm$M, 1, function(xx) {
            mean(y2var_k1[which(xx == 1)])/sum(xx)
        })
        y2var_k2 <- x$y2var_k2
        ay2var_k2 <- apply(rm$M, 1, function(xx) {
            mean(y2var_k2[which(xx == 1)])
        })
    }
    mi.coef <- NULL
    mi.cov <- NULL
    mi.tau2 <- c()
    mi.qe <- c()
    mi.i2 <- c()
    mi.h2 <- c()
    for (j in 1:n.imp) {
        y[unknown] <- x$unknown$y[, j]
        ay <- c(rm.M %*% y)
        if (measure == "smc" || measure == "smd") {
            ay.var <- ay2var_k1 + ay2var_k2 * ay^2
        }
        W_fe <- diag(1/ay.var)
        P_fe <- W_fe - W_fe %*% aX %*% solve(t(aX) %*% W_fe %*% 
            aX) %*% t(aX) %*% W_fe
        tau2_j <- .tau2.reml(tau2, ay, ay.var, aX, maxiter, tol)
        W <- diag(1/(ay.var + tau2_j))
        inv_XtWX <- solve(t(aX) %*% W %*% aX)
        h2_j <- 1 + tau2_j/df * sum(diag(P_fe))
        mi.coef <- cbind(mi.coef, inv_XtWX %*% t(aX) %*% W %*% 
            ay)
        mi.cov <- cbind(mi.cov, c(inv_XtWX))
        mi.tau2 <- c(mi.tau2, tau2_j)
        mi.qe <- c(mi.qe, max(0, t(ay) %*% P_fe %*% ay))
        mi.i2 <- c(mi.i2, max(0, 1 - 1/h2_j))
        mi.h2 <- c(mi.h2, h2_j)
    }
    coef <- apply(mi.coef, 1, mean)
    cov <- .pool.cov(mi.coef, mi.cov)
    se <- sqrt(diag(cov))
    z <- coef/se
    x$model$coef <- coef
    x$model$cov <- cov
    x$model$se <- se
    x$model$z <- z
    x$model$p.value <- 2 * pnorm(-abs(z))
    if (n.coef > 1) {
        coef <- coef[-1]
        cov <- as.matrix(cov[-1, -1])
        qr <- qr(cov)
        pivot <- qr$pivot[1:qr$rank]
        chisq <- c(t(coef[pivot]) %*% solve(cov[pivot, pivot]) %*% 
            coef[pivot])
        df <- length(pivot)
        x$model$q <- data.frame(chisq, df, p.value = 1 - pchisq(chisq, 
            df))
    }
    else {
        x$model$q <- NULL
    }
    f_df2 <- .pool.chi2(mi.qe, df)
    f <- f_df2[1]
    df2 <- f_df2[2]
    x$heterogeneity <- list(tau2 = mean(mi.tau2), h2 = mean(mi.h2), 
        i2 = mean(mi.i2), q = data.frame(f, df1 = df, df2, p.value = 1 - 
            pf(f, df, df2)))
    .linearHypothesis.meta.nsue(x)
}
metalm.meta.nsue <-
function (x, formula, maxiter = 100, tol = 1e-06, ...) 
{
    call <- match.call()
    if (!inherits(x, "meta.nsue")) {
        .stop(call, "It must be a 'meta.nsue' object XXX")
    }
    x$model <- .check.formula(call, formula, x$known$n.stud + 
        x$unknown$n.stud)
    x$hypotheses <- .check.hypotheses(call, NULL, x$model)
    .metalm.meta.nsue(x, maxiter, tol)
}
.meta.nsue <-
function (x, X, model, hypotheses, n.imp, n.bins, maxiter, tol, 
    ...) 
{
    y <- x$y
    measure <- x$measure
    n.stud <- length(y)
    known <- which(!is.na(y))
    known.n.stud <- length(known)
    known.y <- y[known]
    unknown <- which(is.na(y))
    unknown.n.stud <- length(unknown)
    unknown.y_lo <- x$y_lo[unknown]
    unknown.y_up <- x$y_up[unknown]
    if (known.n.stud) {
        initial_coef1 <- mean(known.y)
    }
    else {
        initial_coef1 <- 0
    }
    if (measure == "cor" || measure == "cor in smd") {
        y.var = x$y.var
        known.y.var <- y.var[known]
        unknown.y.var <- y.var[unknown]
    }
    if (measure == "smc" || measure == "smd") {
        y2var_k1 <- x$y2var_k1
        known.y2var_k1 = y2var_k1[known]
        unknown.y2var_k1 = y2var_k1[unknown]
        y2var_k2 <- x$y2var_k2
        known.y2var_k2 = y2var_k2[known]
        unknown.y2var_k2 = y2var_k2[unknown]
        known.y.var <- known.y2var_k1 + known.y2var_k2 * known.y^2
    }
    rm <- x$rm
    rm.M <- t(unname(model.matrix(~0 + x$labels)))
    rm.M <- rm.M[unique(c(1:nrow(rm.M) %*% rm.M)), ]
    rm$n.stud = nrow(rm.M)
    rm$M <- rm.M
    rm$weights <- apply(apply(rm.M, 1, function(x) {
        x/(1 + (sum(x) - 1) * rm$r)
    }), 1, sum)
    known.weights = rm$weights[known]
    unknown.weights = rm$weights[unknown]
    if (measure == "cor" || measure == "cor in smd") {
        mll <- function(tau2_coef, known, known.y, known.y.var, 
            known.weights, unknown, unknown.y_lo, unknown.y_up, 
            unknown.y.var, unknown.weights, X) {
            tau2 <- tau2_coef[1]
            if (tau2 < 0) {
                return(Inf)
            }
            mu <- X %*% tau2_coef[-1]
            known.sigma2 <- known.y.var + tau2
            unknown.mu <- mu[unknown]
            unknown.sigma <- sqrt(unknown.y.var + tau2)
            sum(known.weights * (log(known.sigma2) + (known.y - 
                mu[known])^2/known.sigma2))/2 - sum(unknown.weights * 
                log(pnorm((unknown.y_up - unknown.mu)/unknown.sigma) - 
                  pnorm((unknown.y_lo - unknown.mu)/unknown.sigma)))
        }
    }
    if (measure == "smc" || measure == "smd") {
        mll_d <- function(tau2_coef, known, known.y, known.y.var, 
            known.weights, unknown, unknown.y_lo, unknown.y_lo.var, 
            unknown.y_up, unknown.y_up.var, unknown.weights, 
            X) {
            tau2 <- tau2_coef[1]
            if (tau2 < 0) {
                return(Inf)
            }
            mu <- X %*% tau2_coef[-1]
            known.sigma2 <- known.y.var + tau2
            unknown.mu <- mu[unknown]
            unknown.y_lo.sigma <- sqrt(unknown.y_lo.var + tau2)
            unknown.y_up.sigma <- sqrt(unknown.y_up.var + tau2)
            sum(known.weights * (log(known.sigma2) + (known.y - 
                mu[known])^2/known.sigma2))/2 - sum(unknown.weights * 
                log(pnorm((unknown.y_up - unknown.mu)/unknown.y_up.sigma) - 
                  pnorm((unknown.y_lo - unknown.mu)/unknown.y_lo.sigma)))
        }
    }
    if (measure == "cor" || measure == "cor in smd") {
        tau2_coef <- optim(c(0, initial_coef1, rep(0, ncol(X) - 
            1)), mll, gr = NULL, known, known.y, known.y.var, 
            known.weights, unknown, unknown.y_lo, unknown.y_up, 
            unknown.y.var, unknown.weights, X)$par
    }
    if (measure == "smc" || measure == "smd") {
        tau2_coef <- optim(c(0, initial_coef1, rep(0, ncol(X) - 
            1)), mll_d, gr = NULL, known, known.y, known.y.var, 
            known.weights, unknown, unknown.y_lo, unknown.y2var_k1 + 
                unknown.y2var_k2 * unknown.y_lo^2, unknown.y_up, 
            unknown.y2var_k1 + unknown.y2var_k2 * unknown.y_up^2, 
            unknown.weights, X)$par
    }
    tau2 <- tau2_coef[1]
    mu <- X %*% tau2_coef[-1]
    mi_y <- NULL
    if (unknown.n.stud) {
        if (measure == "cor" || measure == "cor in smd") {
            mi <- function(tau2, mu, y_lo, y_up, y.var, rm.mu, 
                rm.var, rm.y, n.imp) {
                sigma2 <- y.var + tau2
                sigma <- sqrt(sigma2)
                if (length(rm.mu)) {
                  rm.sigma <- sqrt(rm.var + tau2)
                  Sxx <- rm.sigma %*% t(rm.sigma) * (diag(1 - 
                    rm$r, length(rm.mu)) + rm$r)
                  Sxy <- rm.sigma * sigma * rm$r
                  beta <- solve(Sxx) %*% Sxy
                  mus <- rm.y %*% beta
                  sigma <- sqrt(sigma2 - t(beta) %*% Sxx %*% 
                    beta)
                }
                else {
                  mus <- rep(mu, n.imp)
                }
                q <- rep(NA, n.imp)
                to_imp <- 1:n.imp
                while (length(to_imp)) {
                  q[to_imp] <- rnorm(length(to_imp), mus[to_imp], 
                    sigma)
                  to_imp <- which(q <= y_lo & q >= y_up)
                }
                q
            }
        }
        if (measure == "smc" || measure == "smd") {
            mi_d <- function(tau2, mu, y_lo, y_up, y2var_k1, 
                y2var_k2, rm.mu, rm.var_k1, rm.var_k2, rm.y, 
                n.imp, n.bins) {
                sigma2 <- y2var_k1 + y2var_k2 * mu^2 + tau2
                sigma <- sqrt(sigma2)
                if (length(rm.mu)) {
                  rm.sigma <- sqrt(rm.var_k1 + rm.var_k2 * rm.mu^2 + 
                    tau2)
                  Sxx <- rm.sigma %*% t(rm.sigma) * (diag(1 - 
                    rm$r, length(rm.mu)) + rm$r)
                  Sxy <- rm.sigma * sigma * rm$r
                  beta <- solve(Sxx) %*% Sxy
                  mus <- rm.y %*% beta
                  sigma <- sqrt(sigma2 - t(beta) %*% Sxx %*% 
                    beta)
                }
                else {
                  mus <- rep(mu, n.imp)
                }
                width <- (y_up - y_lo)/n.bins
                y <- sample(seq(y_lo + width/2, y_up - width/2, 
                  width))
                q <- c()
                for (imp in 1:n.imp) {
                  raw_dens <- dnorm(y, mus[imp], sigma) * (y2var_k1 + 
                    y2var_k2 * y^2 + tau2)
                  pfun <- cumsum(raw_dens/sum(raw_dens))
                  p <- runif(1)
                  j <- 1
                  while (p > pfun[j]) {
                    j <- j + 1
                  }
                  q <- c(q, y[j])
                }
                q
            }
        }
        for (i in 1:unknown.n.stud) {
            index <- unknown[i]
            rm.indexs <- which(x$labels == x$labels[index] & 
                1:n.stud < index)
            rm.mu = c()
            if (measure == "cor" || measure == "cor in smd") {
                rm.var <- c()
            }
            if (measure == "smc" || measure == "smd") {
                rm.var_k1 <- c()
                rm.var_k2 <- c()
            }
            rm.y = NULL
            for (j in rm.indexs) {
                is.known = !is.na(x$y[j])
                rm.mu <- c(rm.mu, mu[j])
                if (measure == "cor" || measure == "cor in smd") {
                  rm.var <- c(rm.var, y.var[j])
                }
                if (measure == "smc" || measure == "smd") {
                  rm.var_k1 <- c(rm.var_k1, y2var_k1[j])
                  rm.var_k2 <- c(rm.var_k2, y2var_k2[j])
                }
                if (is.known) {
                  rm.y <- cbind(rm.y, rep(x$y[j], n.imp))
                }
                else {
                  rm.y <- cbind(rm.y, mi_y[match(j, unknown), 
                    ])
                }
            }
            if (measure == "cor" || measure == "cor in smd") {
                mi_y <- rbind(mi_y, mi(tau2, mu[index], unknown.y_lo[i], 
                  unknown.y_up[i], unknown.y.var[i], rm.mu, rm.var, 
                  rm.y, n.imp))
            }
            if (measure == "smc" || measure == "smd") {
                mi_y <- rbind(mi_y, mi_d(tau2, mu[index], unknown.y_lo[i], 
                  unknown.y_up[i], y2var_k1[index], y2var_k2[index], 
                  rm.mu, rm.var_k1, rm.var_k2, rm.y, n.imp, n.bins))
            }
        }
        colnames(mi_y) <- NULL
    }
    else {
        n.imp <- 0
    }
    obj <- list(measure = measure, known = list(n.stud = length(known), 
        i = known, y = known.y), unknown = list(n.stud = unknown.n.stud, 
        i = unknown, n.imp = n.imp, y = mi_y))
    if (measure == "cor" || measure == "cor in smd") {
        obj$y.var <- y.var
    }
    if (measure == "smc" || measure == "smd") {
        obj$y2var_k1 <- y2var_k1
        obj$y2var_k2 <- y2var_k2
    }
    if (measure == "cor in smd") {
        obj$smd <- x$smd
    }
    obj$labels = x$labels
    obj$rm = rm
    obj$heterogeneity = list(tau2 = tau2)
    obj$model = model
    obj$hypotheses = hypotheses
    class(obj) <- "meta.nsue"
    .metalm.meta.nsue(obj, maxiter, tol)
}
meta.nsue <-
function (x, data = data.frame(), formula = ~1, hypotheses = NULL, 
    n.imp = 50, n.bins = 200, maxiter = 200, tol = 1e-06, ...) 
{
    call <- match.call()
    if (!inherits(x, "nsue")) {
        .stop(call, "Use an smc_from_t, smd_from_t or r_from_z call as the first (nsue) argument.")
    }
    n.stud <- length(x$y)
    X <- .check.data(call, data, n.stud)
    model <- .check.formula(call, formula, n.stud)
    hypotheses <- .check.hypotheses(call, hypotheses, model)
    if (n.imp < 2) {
        .stop(call, "At least 2 imputations XXX")
    }
    .meta.nsue(x, X, model, hypotheses, n.imp, n.bins, maxiter, 
        tol)
}
plot.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    forest.meta.nsue(x)
}
.pool.chi2 <-
function (chi2, df1) 
{
    m <- length(chi2)
    r <- (1 + 1/m) * var(sqrt(chi2))
    c((mean(chi2)/df1 - (m + 1)/(m - 1) * r)/(1 + r), (m - 1)/df1^(3/m) * 
        (1 + 1/r)^2)
}
.pool.cov <-
function (x, x_cov) 
{
    cov0 <- matrix(apply(x_cov, 1, mean), nrow(x))
    var0 <- diag(cov0)
    var.increase <- sqrt((var0 + (1 + 1/ncol(x)) * apply(x, 1, 
        var))/var0)
    var.increase %*% t(var.increase) * cov0
}
.pool.var <-
function (x, x_var) 
{
    if (is.vector(x)) {
        return(mean(x_var) + (1 + 1/length(x)) * var(x))
    }
    apply(x_var, 1, mean) + (1 + 1/ncol(x)) * apply(x, 1, var)
}
.print.description <-
function (x) 
{
    cat("\n")
    cat("Meta-analysis description:\n")
    cat("- Measure:", .format.measure(x$measure), "\n")
    cat("- Known measures:", x$known$n.stud, "\n")
    if (x$known$n.stud == 0) {
        .warning("No known measures!")
    }
    cat("- Non-statistically-significant unknown measures:", 
        x$unknown$n.stud, "\n")
    if (x$rm$n.stud < x$known$n.stud + x$unknown$n.stud) {
        cat("- Measures after combination of repeated-measures:", 
            x$rm$n.stud, "\n")
    }
    cat("- Imputations:", x$unknown$n.imp, "\n")
    cat("- Model: measure", x$model$formula, "\n")
    cat("\n")
}
.print.heterogeneity_and_model <-
function (x) 
{
    measure <- x$measure
    cat("Residual heterogeneity (tau^2):", .format.4pos(x$heterogeneity$tau2), 
        "  I^2:", paste(.format.perc2(x$heterogeneity$i2), "%", 
            sep = ""), "  H^2:", .format.2pos(x$heterogeneity$h2), 
        "\n")
    cat("F-statistic (heterogeneity):", .format.2pos(x$heterogeneity$q$f), 
        "on", x$heterogeneity$q$df1, "and", .format.1pos(x$heterogeneity$q$df2), 
        "df  Pr(>F):", .format.prob(x$heterogeneity$q$p.value), 
        "\n")
    cat("\n")
    cat("Model:\n")
    coef <- x$model$coef
    se <- x$model$se
    p.value <- x$model$p.value
    ci.low <- coef + qnorm(0.025) * se
    ci.up <- coef + qnorm(0.975) * se
    sign <- .format.sign(p.value)
    if (measure == "cor" || measure == "cor in smd") {
        table <- cbind(cbind(.format.3(tanh(coef)), .format.4pos(se), 
            .format.4(x$model$z), .format.prob(p.value), .format.3(tanh(cbind(ci.low, 
                ci.up)))), sign)
        colnames(table) <- c("Corr", "Std. Error", "z value", 
            "Pr(>|z|)", "CI(low)", "CI(up)", "")
    }
    if (measure == "smc" || measure == "smd") {
        table <- cbind(cbind(.format.4(coef), .format.4pos(se), 
            .format.4(x$model$z), .format.prob(p.value), .format.4(cbind(ci.low, 
                ci.up))), sign)
        colnames(table) <- c("Estimate", "Std. Error", "z value", 
            "Pr(>|z|)", "CI(low)", "CI(up)", "")
    }
    rownames(table) <- substr(x$model$labels, 1, 12)
    print(table, quote = FALSE, right = TRUE, print.gap = 2)
    if (is.null(x$model$q)) {
        cat("Chisq-statistic (model): NA\n")
    }
    else {
        cat("Chisq-statistic (model):", .format.2pos(x$model$q$chisq), 
            "on", x$model$q$df, "df  Pr(>chisq):", .format.prob(x$model$q$p.value), 
            "\n")
    }
    cat("\n")
}
print.leave1out.nsue <-
function (x, ...) 
{
    if (!inherits(x, "leave1out.nsue")) {
        .stop(match.call(), "The argument must be a 'leave1out.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis description:\n")
    cat("- Measure:", .format.measure(x[[1]]$meta.nsue$measure), 
        "\n")
    cat("- Model: measure", x[[1]]$meta.nsue$model$formula, "\n")
    cat("\n")
    for (i in 1:length(x)) {
        cat("\n")
        cat("Discarded study:", x[[i]]$study, "\n")
        cat("\n")
        .print.heterogeneity_and_model(x[[i]]$meta.nsue)
        .summary.meta.nsue(x[[i]]$meta.nsue)
    }
    .print.sign()
    cat("\n")
}
print.meta.nsue <-
function (x, ...) 
{
    if (!inherits(x, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis description:\n")
    cat("- Measure:", .format.measure(x$measure), "\n")
    cat("- Known measures:", x$known$n.stud, "\n")
    if (x$known$n.stud == 0) {
        .warning("No known measures!")
    }
    cat("- Non-statistically-significant unknown measures:", 
        x$unknown$n.stud, "\n")
    if (x$rm$n.stud < x$known$n.stud + x$unknown$n.stud) {
        cat("- Measures after combination of repeated-measures:", 
            x$rm$n.stud, "\n")
    }
    cat("- Imputations:", x$unknown$n.imp, "\n")
    cat("- Model: measure", x$model$formula, "\n")
    cat("\n")
    .print.heterogeneity_and_model(x)
    .summary.meta.nsue(x)
    .print.sign()
    cat("\n")
}
print.nsue <-
function (x, ...) 
{
    cat("\n")
    cat("'nsue' object description:\n")
    cat("- Measure:", .format.measure(x$measure), "\n")
    known.n.stud = sum(!is.na(x$y))
    unknown.n.stud = sum(is.na(x$y))
    cat("- Known measures:", known.n.stud, "\n")
    if (known.n.stud == 0) {
        .warning("No known measures!")
    }
    cat("- Non-statistically-significant unknown measures:", 
        sum(is.na(x$y)), "\n")
    cat("\n")
}
.print.sign <-
function () 
{
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}
residuals.meta.nsue <-
function (object, ...) 
{
    fitted <- fitted(object)
    known.i <- object$known$i
    residuals <- c()
    residuals[known.i] <- object$known$y - fitted[known.i]
    unknown.i <- object$unknown$i
    if (object$unknown$n.stud) {
        residuals[unknown.i] <- apply(object$unknown$y, 1, mean) - 
            fitted[unknown.i]
    }
    residuals
}
.signif.up <-
function (x, digits = 6) 
{
    power <- 10^(round(digits) - 1 - floor(log10(abs(x))))
    ceiling(x * power)/power
}
smc_from_t <-
function (t, n, alpha = 0.05, labels = "study", rm.r = 0.3) 
{
    call <- match.call()
    if (missing(t) || missing(n)) {
        .stop(call, "You must specify t and n")
    }
    if (!is.numeric(t)) {
        .stop(call, "t is not a numeric vector")
    }
    n.stud <- length(t)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n <- .check.n(call, n, 3, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    labels <- .check.labels(call, labels, n.stud)
    for (i in 1:n.stud) {
        if ((is.na(t[i]) && is.na(alpha[i])) || is.na(n[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    df <- n - 1
    inv_n <- 1/n
    j <- .d_j(df)
    k_t2d <- j * sqrt(inv_n)
    y_up <- k_t2d * qt(1 - alpha/2, df)
    obj <- list(measure = "smc", y = k_t2d * t, y_lo = -y_up, 
        y_up = y_up, y2var_k1 = inv_n, y2var_k2 = 1 - (df - 2)/(df * 
            j^2), labels = labels, rm = list(r = rm.r))
    class(obj) <- "nsue"
    obj
}
smd_from_t <-
function (t, n1, n2, alpha = 0.05, labels = "study", rm.r = 0.3) 
{
    call <- match.call()
    if (missing(t) || missing(n1) || missing(n2)) {
        .stop(call, "You must specify t, n1 and n2")
    }
    if (!is.numeric(t)) {
        .stop(call, "t is not a numeric vector")
    }
    n.stud <- length(t)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n1 <- .check.n(call, n1, 2, n.stud)
    n2 <- .check.n(call, n2, 2, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    labels <- .check.labels(call, labels, n.stud)
    for (i in 1:n.stud) {
        if ((is.na(t[i]) && is.na(alpha[i])) || is.na(n1[i]) || 
            is.na(n2[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    n <- n1 + n2
    df <- n - 2
    inv_n1_n2 <- 1/n1 + 1/n2
    j <- .d_j(df)
    k_t2d <- j * sqrt(inv_n1_n2)
    y_up = k_t2d * qt(1 - alpha/2, df)
    obj <- list(measure = "smd", y = k_t2d * t, y_lo = -y_up, 
        y_up = y_up, y2var_k1 = inv_n1_n2, y2var_k2 = 1 - (df - 
            2)/(df * j^2), labels = labels, rm = list(r = rm.r))
    class(obj) <- "nsue"
    obj
}
.stop <-
function (call, message) 
{
    cat("\n")
    print(call)
    stop(paste0(message, "\n "), call. = FALSE)
}
summary.leave1out.nsue <-
function (object, ...) 
{
    if (!inherits(object, "leave1out.nsue")) {
        .stop(match.call(), "The argument must be a 'leave1out.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis model:", object[[1]]$meta.nsue$measure, 
        object[[1]]$meta.nsue$model$formula, "\n")
    for (i in 1:length(object)) {
        cat("\n")
        cat("Discarded study:", object[[i]]$study, "\n")
        cat("\n")
        .summary.meta.nsue(object[[i]]$meta.nsue)
    }
    .print.sign()
    cat("\n")
    invisible(object)
}
.summary.meta.nsue <-
function (x) 
{
    one <- which(!is.na(x$hypotheses$z))
    if (length(one)) {
        measure <- x$measure
        cat("One-row hypotheses:\n")
        coef <- x$hypotheses$coef[one]
        se <- x$hypotheses$se[one]
        p.value <- x$hypotheses$p.value[one]
        ci.low <- coef + qnorm(0.025) * se
        ci.up <- coef + qnorm(0.975) * se
        sign <- .format.sign(p.value)
        if (measure == "cor" || measure == "cor in smd") {
            hypotheses <- cbind(cbind(.format.3(tanh(coef)), 
                .format.4pos(se), .format.4(x$hypotheses$z[one]), 
                .format.prob(p.value), .format.3(tanh(cbind(ci.low, 
                  ci.up)))), sign)
            colnames(hypotheses) <- c("Corr", "Std. Error", "z value", 
                "Pr(>|z|)", "CI(low)", "CI(up)", "")
        }
        if (measure == "smc" || measure == "smd") {
            hypotheses <- cbind(cbind(.format.4(coef), .format.4pos(se), 
                .format.4(x$hypotheses$z[one]), .format.prob(p.value), 
                .format.4(cbind(ci.low, ci.up))), sign)
            colnames(hypotheses) <- c("Estimate", "Std. Error", 
                "z value", "Pr(>|z|)", "CI(low)", "CI(up)", "")
        }
        rownames(hypotheses) <- substr(x$hypotheses$labels[one], 
            1, 12)
        print(hypotheses, quote = FALSE, right = TRUE, print.gap = 2)
        cat("\n")
    }
    multi <- which(!is.na(x$hypotheses$chisq))
    if (length(multi)) {
        cat("Multi-row hypotheses:\n")
        p.value <- x$hypotheses$p.value[multi]
        hypotheses <- cbind(.format.2pos(x$hypotheses$chisq[multi]), 
            .format.0pos(x$hypotheses$df[multi]), .format.prob(p.value), 
            .format.sign(p.value))
        colnames(hypotheses) <- c("chisq", "df", "Pr(>chisq)", 
            "")
        rownames(hypotheses) <- substr(x$hypotheses$labels[multi], 
            1, 12)
        print(hypotheses, quote = FALSE, right = TRUE, print.gap = 2)
        cat("\n")
    }
}
summary.meta.nsue <-
function (object, ...) 
{
    if (!inherits(object, "meta.nsue")) {
        .stop(match.call(), "The argument must be a 'meta.nsue' object")
    }
    cat("\n")
    cat("Meta-analysis model:", object$measure, object$model$formula, 
        "\n")
    cat("\n")
    .summary.meta.nsue(object)
    .print.sign()
    cat("\n")
    invisible(object)
}
.tau2.reml <-
function (tau2, y, y.var, X, maxiter, tol) 
{
    for (iter in 1:maxiter) {
        old_tau2 <- tau2
        W <- diag(1/(y.var + tau2))
        P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% 
            W
        tau2 <- max(0, tau2 + solve(sum(diag(P %*% P))) %*% (t(y) %*% 
            P %*% P %*% y - sum(diag((P)))))
        if (abs(tau2 - old_tau2) < tol) {
            break
        }
    }
    tau2
}
.warning <-
function (message) 
{
    warning(message, call. = FALSE)
}
z_from_r <-
function (r, n, alpha = 0.05, labels = "study", rm.r = 0.3) 
{
    call <- match.call()
    if (missing(r) || missing(n)) {
        .stop(call, "You must specify r and n")
    }
    if (!is.numeric(r)) {
        .stop(call, "r is not a numeric vector")
    }
    if (any(r < -1, na.rm = TRUE) || any(r > 1, na.rm = TRUE)) {
        .stop(call, "r cannot be <= -1 or >= 1")
    }
    n.stud <- length(r)
    if (!n.stud) {
        .stop(call, "No studies to meta-analyze")
    }
    n <- .check.n(call, n, 4, n.stud)
    alpha <- .check.alpha(call, alpha, n.stud)
    labels <- .check.labels(call, labels, n.stud)
    for (i in 1:n.stud) {
        if ((is.na(r[i]) && is.na(alpha[i])) || is.na(n[i])) {
            stop("Not enough information in study", labels[i])
        }
    }
    y_up = atanh((1 + (n - 2)/qt(alpha/2, n - 2)^2)^-0.5)
    obj <- list(measure = "cor", y = atanh(r), y_lo = -y_up, 
        y_up = y_up, y.var = 1/(n - 3), labels = labels, rm = list(r = rm.r))
    class(obj) <- "nsue"
    obj
}
