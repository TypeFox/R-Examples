expectreg.qp <-
function (formula, data = NULL, id = NA, smooth = c("schall", 
    "acv", "fixed"), lambda = 1, expectiles = NA) 
{
    smooth = match.arg(smooth)
    if (any(is.na(expectiles)) || !is.vector(expectiles) || any(expectiles > 
        1) || any(expectiles < 0)) {
        p <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
            0.98, 0.99)
        row.grid = 3
        col.grid = 4
    }
    else {
        p <- expectiles
        row.grid = floor(sqrt(length(p)))
        col.grid = ceiling(sqrt(length(p)))
        if (length(p) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    y = eval(as.expression(formula[[2]]), envir = data, enclos = environment(formula))
    attr(y, "name") = deparse(formula[[2]])
    n <- length(y)
    mp <- length(p)
    i.ja <- 0
    if (attr(terms(formula), "intercept") == "1") 
        i.ja <- 1
    types = list()
    design = list()
    helper = list()
    x = list()
    if (attr(terms(formula), "intercept") == "1") {
        design[[1]] = rb(matrix(1, nrow = n, ncol = 1), "parametric", 
            center = F)
        types[[1]] = "intercept"
    }
    else if (formula[[3]] == ".") {
        design[[1]] = rb(data[, names(data) != all.vars(formula[[2]])], 
            "parametric")
    }
    if (length(labels(terms(formula))) > 0) 
        for (i in 1:length(labels(terms(formula)))) {
            types[[i.ja + i]] = strsplit(strsplit(labels(terms(formula))[i], 
                "(", fixed = TRUE)[[1]][2], ",", fixed = TRUE)[[1]][2]
            if (is.na(types[[i.ja + i]])) 
                types[[i.ja + i]] = "pspline"
            if (types[[i.ja + i]] == "pspline" || types[[i.ja + 
                i]] == " \"pspline\")") 
                types[[i.ja + i]] <- "pspline"
            if (types[[i.ja + i]] == " \"parametric\")" || types[[i.ja + 
                i]] == "parametric") {
                lambda[i.ja + i] <- 0
                P <- diag(0, 2, 2)
                types[[i.ja + i]] <- "parametric"
                x.mom <- eval(parse(text = labels(terms(formula))[i], 
                  srcfile = NULL), envir = data, enclos = environment(formula))
                if (is.factor(x.mom$x)) 
                  design[[i.ja + i]] = x.mom
                else {
                  x.mom = x.mom$x
                  knots <- c(min(x.mom) - 2, min(x.mom) - 1e-04, 
                    max(x.mom) + 1e-04, max(x.mom) + 2)
                  B <- splineDesign(knots = knots, x = x.mom, 
                    ord = 2)
                  design[[i.ja + i]] <- list(B = B, P = P, x = x.mom, 
                    type = types[[i.ja + i]], constraint = matrix(0, 
                      nrow = 2, ncol = ncol(B)))
                  design[[i.ja + i]]$xname <- eval(parse(text = labels(terms(formula))[i], 
                    srcfile = NULL), envir = data, enclos = environment(formula))$xname
                }
            }
            else if (strsplit(labels(terms(formula))[i], "(", 
                fixed = TRUE)[[1]][1] != "mono") {
                if (sub("center=TRUE", "center=FALSE", labels(terms(formula))[i]) == 
                  labels(terms(formula))[i]) {
                  basetext = paste("rb(center=FALSE,", strsplit(labels(terms(formula))[i], 
                    "rb(", fixed = TRUE)[[1]][2], sep = "")
                }
                else {
                  basetext = sub("center=TRUE", "center=FALSE", 
                    labels(terms(formula))[i])
                }
                design[[i.ja + i]] = eval(parse(text = basetext), 
                  envir = data, enclos = environment(formula))
            }
            else design[[i.ja + i]] = eval(parse(text = labels(terms(formula))[i]), 
                envir = data, enclos = environment(formula))
        }
    Bx = design[[1]][[1]]
    K <- vector(length = length(design))
    K[1] <- ncol(Bx)
    P = design[[1]][[2]]
    P = t(P) %*% P
    P.list <- list()
    P.list[[1]] <- P
    bigP <- P
    x[[1]] = design[[1]][[3]]
    names(x)[1] = design[[1]]$xname[1]
    if (types[[1]] == "markov") 
        helper[[1]] = list(design[[1]][[5]], design[[1]][[6]])
    else if (types[[1]] == "krig") 
        helper[[1]] = list(design[[1]][[7]])
    else helper[[1]] = list(NA)
    Cmat = list()
    Cmat[[1]] = as.matrix(design[[1]]$constraint)
    for (i in 2:length(p)) {
        Cmat[[1]] = rbind(cbind(Cmat[[1]], matrix(0, nrow = nrow(Cmat[[1]]), 
            ncol = ncol(as.matrix(design[[1]]$constraint)))), 
            cbind(matrix(0, nrow = nrow(as.matrix(design[[1]]$constraint)), 
                ncol = ncol(Cmat[[1]])), as.matrix(design[[1]]$constraint)))
    }
    if (length(design) > 1) {
        for (i in 2:length(design)) {
            Bx = cbind(Bx, design[[i]][[1]])
            K[i] <- ncol(design[[i]][[1]])
            types[[i]] = design[[i]][[4]]
            bigP = rbind(cbind(P, matrix(0, nrow = nrow(P), ncol = ncol(design[[i]][[2]]))), 
                cbind(matrix(0, nrow = nrow(design[[i]][[2]]), 
                  ncol = ncol(P)), design[[i]][[2]]))
            P.list[[i]] <- t(design[[i]][[2]]) %*% design[[i]][[2]]
            x[[i]] = design[[i]][[3]]
            names(x)[i] = design[[i]]$xname[1]
            if (types[[i]] == "markov") 
                helper[[i]] = list(design[[i]][[5]], design[[i]][[6]])
            else if (types[[i]] == "krig") 
                helper[[i]] = list(design[[i]][[7]])
            else helper[[i]] = list(NA)
            Cmat[[i]] = as.matrix(design[[i]]$constraint)
            for (k in 2:length(p)) {
                Cmat[[i]] = rbind(cbind(Cmat[[i]], matrix(0, 
                  nrow = nrow(Cmat[[i]]), ncol = ncol(as.matrix(design[[i]]$constraint)))), 
                  cbind(matrix(0, nrow = nrow(as.matrix(design[[i]]$constraint)), 
                    ncol = ncol(Cmat[[i]])), as.matrix(design[[i]]$constraint)))
            }
        }
    }
    if (all(is.na(id))) 
        nc.bf = qpvzentr(design, lambda, types, p, y, P.list, 
            Cmat, Bx, K, i.ja, smooth)
    else {
        id = eval(id, envir = data, enclos = environment(formula))
        nc.bf = qplongi(design, lambda, types, p, y, P.list, 
            Cmat, Bx, K, i.ja, smooth, id)
    }
    nc.bf$fitted <- nc.bf$fitted + nc.bf$alpha
    add_alpha <- function(x) x <- x + nc.bf$alpha
    Z <- lapply(nc.bf$values, add_alpha)
    alpha <- nc.bf$alpha
    if (attr(terms(formula), "intercept") == "1") {
        nc.bf$fitted[, -1] <- nc.bf$fitted[, -1] + nc.bf$fitted[, 
            1]
        add_intercept <- function(x) x <- x + nc.bf$fitted[, 
            1]
        Z <- lapply(nc.bf$values, add_intercept)
        Z[[1]] <- nc.bf$values[[1]]
    }
    coefficients <- nc.bf$coefficients
    coeff.vec <- nc.bf$coeff.vec
    dfi <- nc.bf$dfi
    final.lambdas <- nc.bf$lambdanp
    final.lambdas = list(final.lambdas)
    if (i.ja == 1) {
        intercepts <- Z[[1]][1, ] + alpha
        Z <- Z[-1]
        types <- types[-1]
        x <- x[-1]
        coefficients <- coefficients[-1]
        design <- design[-1]
    }
    else intercepts <- rep(0, length(p))
    for (j in 1:length(design)) {
        names(Z)[j] = design[[j]]$xname[1]
        names(coefficients)[j] = design[[j]]$xname[1]
    }
    yyhat = matrix(nc.bf$fitted, nrow = n, ncol = length(p))
    ncexpect <- list(lambda = final.lambdas, intercepts = intercepts, 
        values = Z, coefficients = coefficients, response = y, 
        formula = formula, asymmetries = p, effects = types, 
        helper = helper, covariates = x, design = Bx, bases = design, 
        fitted = yyhat, weights = nc.bf$weights, part.resid = nc.bf$part.resid, 
        dfi = dfi, coeff.vec = coeff.vec, alpha = alpha, tau2 = nc.bf$tau2, 
        sig2 = nc.bf$sig2)
    ncexpect$predict <- function(newdata = NULL) {
        BB = list()
        values = list()
        bmat = NULL
        betavec = NULL
        fitted <- matrix(NA, ncol = length(p), nrow = length(newdata))
        it = 1
        for (k in it:length(coefficients)) {
            BB[[k]] = predict(design[[k]], newdata)
            values[[k]] <- BB[[k]] %*% coefficients[[k]]
            bmat = cbind(bmat, BB[[k]])
            betavec = rbind(betavec, coefficients[[k]])
        }
        add_alpha <- function(x) x <- x + nc.bf$alpha
        values <- lapply(values, add_alpha)
        names(values) = names(coefficients)
        if (attr(terms(formula), "intercept") == "1") {
            bmat = cbind(1, bmat)
            betavec = rbind(intercepts, betavec)
        }
        fitted = bmat %*% betavec + alpha
        names(values) = names(coefficients)
        list(fitted = fitted, values = values)
    }
    class(ncexpect) = c("expectreg", "noncross")
    ncexpect
}
