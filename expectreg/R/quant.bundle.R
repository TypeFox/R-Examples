quant.bundle <-
function (formula, data = NULL, smooth = c("schall", "acv", "fixed"), 
    lambda = 1, quantiles = NA, simple = TRUE) 
{
    smooth = match.arg(smooth)
    if (!is.na(charmatch(quantiles[1], "density")) && charmatch(quantiles[1], 
        "density") > 0) {
        pp <- seq(0.01, 0.99, by = 0.01)
    }
    else if (any(is.na(quantiles)) || !is.vector(quantiles) || 
        any(quantiles > 1) || any(quantiles < 0)) {
        pp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 
            0.98, 0.99)
        pp.plot <- 1:length(pp)
        row.grid = 3
        col.grid = 4
    }
    else {
        pp <- quantiles
        pp.plot <- 1:length(pp)
        row.grid = floor(sqrt(length(pp)))
        col.grid = ceiling(sqrt(length(pp)))
        if (length(pp) > row.grid * col.grid) 
            row.grid = row.grid + 1
    }
    np <- length(pp)
    yy = eval(parse(text = formula[2]), envir = data, enclos = environment(formula))
    attr(yy, "name") = deparse(formula[[2]])
    m = length(yy)
    design = list()
    x = list()
    types = list()
    bnd = list()
    Zspathelp = list()
    nb = vector()
    krig.phi = list()
    center = TRUE
    varying = list()
    if (formula[[3]] == "1") {
        design[[1]] = rb(matrix(1, nrow = m, ncol = 1), "parametric", 
            center = F)
        smooth = "fixed"
    }
    else if (formula[[3]] == ".") {
        design[[1]] = rb(data[, names(data) != all.vars(formula[[2]])], 
            "parametric")
        smooth = "fixed"
    }
    else for (i in 1:length(labels(terms(formula)))) {
        types[[i]] = strsplit(labels(terms(formula))[i], "(", 
            fixed = TRUE)[[1]][1]
        if (types[[i]] == labels(terms(formula))[i]) {
            design[[i]] = rb(matrix(eval(parse(text = labels(terms(formula))[i]), 
                envir = data, enclos = environment(formula)), 
                nrow = m), "parametric")
            types[[i]] = "parametric"
            design[[i]]$xname = labels(terms(formula))[i]
        }
        else design[[i]] = eval(parse(text = labels(terms(formula))[i]), 
            envir = data, enclos = environment(formula))
    }
    nterms = length(design)
    varying[[1]] = design[[1]][[9]]
    if (any(!is.na(varying[[1]]))) 
        B = design[[1]][[1]] * varying[[1]]
    else B = design[[1]][[1]]
    DD = as.matrix(design[[1]][[2]])
    x[[1]] = design[[1]][[3]]
    names(x)[1] = design[[1]]$xname
    types[[1]] = design[[1]][[4]]
    bnd[[1]] = design[[1]][[5]]
    Zspathelp[[1]] = design[[1]][[6]]
    nb[1] = ncol(design[[1]][[1]])
    krig.phi[[1]] = design[[1]][[7]]
    center = center && design[[1]][[8]]
    if (length(design) > 1) 
        for (i in 2:length(labels(terms(formula)))) {
            varying[[i]] = design[[i]][[9]]
            if (any(!is.na(varying[[i]]))) 
                B = cbind(B, design[[i]][[1]] * varying[[i]])
            else B = cbind(B, design[[i]][[1]])
            design[[i]][[2]] = as.matrix(design[[i]][[2]])
            DD = rbind(cbind(DD, matrix(0, nrow = nrow(DD), ncol = ncol(design[[i]][[2]]))), 
                cbind(matrix(0, nrow = nrow(design[[i]][[2]]), 
                  ncol = ncol(DD)), design[[i]][[2]]))
            x[[i]] = design[[i]][[3]]
            names(x)[i] = design[[i]]$xname
            types[[i]] = design[[i]][[4]]
            bnd[[i]] = design[[i]][[5]]
            Zspathelp[[i]] = design[[i]][[6]]
            nb[i] = ncol(design[[i]][[1]])
            krig.phi[[i]] = design[[i]][[7]]
            center = center && design[[i]][[8]]
        }
    if (center) {
        B = cbind(1, B)
        DD = rbind(0, cbind(0, DD))
    }
    if (simple) 
        ex = expectreg.ls(formula = formula, data = data, smooth = smooth, 
            lambda = lambda, expectiles = "density", estimate = "restricted")
    else ex = expectreg.ls(formula = formula, data = data, smooth = smooth, 
        lambda = lambda, expectiles = "density", estimate = "bundle")
    bunden = cdf.bundle(ex)
    lala <- matrix(lambda, nrow = nterms, ncol = 2, dimnames = list(1:nterms, 
        c("mean", "residual")))
    vector.a.ma.schall <- matrix(NA, nrow = sum(nb) + (1 * center), 
        ncol = np)
    wq <- rep(NA, length(pp))
    g = bunden$density * (bunden$x[2] - bunden$x[1])
    for (i in 1:length(pp)) {
        wq[i] <- which.min(abs(cumsum(g) - pp[i]))
    }
    hq <- bunden$x[wq]
    ex$asymmetry = hq
    for (q in 1:np) {
        vector.a.ma.schall[, q] = ex$trend.coef + ex$residual.coef * 
            hq[q]
    }
    Z <- list()
    coefficients <- list()
    final.lambdas <- list()
    helper <- list()
    ex$fitted = B %*% vector.a.ma.schall
    if (center) {
        intercept = vector.a.ma.schall[1, ]
        B = B[, -1, drop = FALSE]
        vector.a.ma.schall = vector.a.ma.schall[-1, , drop = FALSE]
    }
    else intercept = rep(0, np)
    for (k in 1:length(design)) {
        final.lambdas[[k]] = lala[k, ]
        names(final.lambdas)[k] = design[[k]]$xname
        partbasis = (sum(nb[0:(k - 1)]) + 1):(sum(nb[0:k]))
        if (types[[k]] == "pspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "markov") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = list(bnd[[k]], Zspathelp[[k]])
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "2dspline") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "radial") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "krig") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = krig.phi[[k]]
            for (i in 1:np) {
                Z[[k]][, i] = design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "random") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "ridge") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "parametric") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        else if (types[[k]] == "special") {
            Z[[k]] <- matrix(NA, m, np)
            coefficients[[k]] = matrix(NA, nrow = nb[k], ncol = np)
            helper[[k]] = NA
            for (i in 1:np) {
                Z[[k]][, i] <- design[[k]][[1]] %*% vector.a.ma.schall[partbasis, 
                  i, drop = FALSE] + intercept[i]
                coefficients[[k]][, i] = vector.a.ma.schall[partbasis, 
                  i, drop = FALSE]
            }
        }
        names(Z)[k] = design[[k]]$xname
        names(coefficients)[k] = design[[k]]$xname
    }
    ex$coefficients = coefficients
    ex$values = Z
    ex$asymmetries = pp
    ex$intercepts = intercept
    ex$predict <- function(newdata = NULL) {
        BB = list()
        values = list()
        bmat = NULL
        for (k in 1:length(coefficients)) {
            BB[[k]] = predict(design[[k]], newdata)
            values[[k]] <- BB[[k]] %*% coefficients[[k]]
            values[[k]] = t(apply(values[[k]], 1, function(x) {
                x + intercept
            }))
            bmat = cbind(bmat, BB[[k]])
        }
        if (center) 
            bmat = cbind(1, bmat)
        fitted = bmat %*% rbind(intercept, vector.a.ma.schall)
        names(values) = names(coefficients)
        list(fitted = fitted, values = values)
    }
    ex
}
