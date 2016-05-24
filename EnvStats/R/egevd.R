egevd <-
function (x, method = "mle", pwme.method = "unbiased", tsoe.method = "med", 
    plot.pos.cons = c(a = 0.35, b = 0), ci = FALSE, ci.parameter = "location", 
    ci.type = "two-sided", ci.method = "normal.approx", information = "observed", 
    conf.level = 0.95) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 3 || length(unique(x)) < 3) 
        stop("'x' must contain at least 3 non-missing distinct values.")
    method <- match.arg(method, c("mle", "pwme", "tsoe"))
    ret.list.method <- method
    pwme.method <- match.arg(pwme.method, c("unbiased", "plotting.position"))
    tsoe.method <- match.arg(tsoe.method, c("med", "lms", "lts"))
    ci.parameter <- match.arg(ci.parameter, c("location", "scale", 
        "shape"))
    information <- match.arg(information, c("observed", "expected"))
    if (method != "tsoe") {
        b0 <- mean(x)
        b1 <- pwMoment(x = x, j = 1, method = pwme.method, plot.pos.cons = plot.pos.cons)
        b2 <- pwMoment(x = x, j = 2, method = pwme.method, plot.pos.cons = plot.pos.cons)
        con <- (2 * b1 - b0)/(3 * b2 - b0) - log(2)/log(3)
        k.init <- 7.859 * con + 2.9554 * con^2
        fcn.to.min <- function(k, b0, b1, b2) {
            ((3 * b2 - b0)/(2 * b1 - b0) - (1 - 3^(-k))/(1 - 
                2^(-k)))^2
        }
        shape.pwme <- nlminb(start = k.init, objective = fcn.to.min, 
            b0 = b0, b1 = b1, b2 = b2)$par
        scale.pwme <- ((2 * b1 - b0) * shape.pwme)/(gamma(1 + 
            shape.pwme) * (1 - 2^(-shape.pwme)))
        location.pwme <- b0 + (scale.pwme * (gamma(1 + shape.pwme) - 
            1))/shape.pwme
    }
    switch(method, pwme = {
        dist.params <- c(location = location.pwme, scale = scale.pwme, 
            shape = shape.pwme)
        if (pwme.method == "unbiased") ret.list.method <- paste("Unbiased", 
            method) else ret.list.method <- paste("Plotting-position ", 
            method, "\n", space(33), "with coefficients\n", space(33), 
            "(a, b) = (", plot.pos.cons["a"], ", ", plot.pos.cons["b"], 
            ")", sep = "")
        if (ci) {
            alpha <- scale.pwme
            k <- shape.pwme
            G <- function(x, k) {
                gaussian.hypergeometric(k, 2 * k, 1 + k, -x)
            }
            v.rr <- function(r, alpha, k) {
                (alpha/(k * (r + 1)^k))^2 * (gamma(1 + 2 * k) * 
                  G(r/(r + 1), k) - gamma(1 + k)^2)
            }
            v.r.rplus1 <- function(r, alpha, k) {
                0.5 * (alpha/k)^2 * ((r + 2)^(-2 * k) * gamma(1 + 
                  2 * k) * G(r/(r + 2), k) + (r + 1)^(-k) * ((r + 
                  1)^(-k) - 2 * (r + 2)^(-k)) * gamma(1 + k)^2)
            }
            v.r.rpluss <- function(r, s, alpha, k) {
                0.5 * (alpha/k)^2 * ((r + s + 1)^(-2 * k) * gamma(1 + 
                  2 * k) * G(r/(r + s + 1), k) - (r + s)^(-2 * 
                  k) * gamma(1 + 2 * k) * G((r + 1)/(r + s), 
                  k) + 2 * (r + 1)^(-k) * ((r + s)^(-k) - (r + 
                  s + 1)^(-k)) * gamma(1 + k)^2)
            }
            v00 <- v.rr(r = 0, alpha = alpha, k = k)
            v01 <- v.r.rplus1(r = 0, alpha = alpha, k = k)
            v02 <- v.r.rpluss(r = 0, s = 2, alpha = alpha, k = k)
            v11 <- v.rr(r = 1, alpha = alpha, k = k)
            v12 <- v.r.rplus1(r = 1, alpha = alpha, k = k)
            v22 <- v.rr(r = 2, alpha = alpha, k = k)
            V <- matrix(c(v00, v01, v02, v01, v11, v12, v02, 
                v12, v22), nrow = 3, byrow = TRUE)
            G.inv <- matrix(0, nrow = 3, ncol = 3)
            v1 <- 1:3
            v2 <- v1^(-k)
            c1 <- gamma(1 + k)
            c2 <- derivative(gamma, 1 + k)
            G.inv[, 1] <- 1/v1
            G.inv[, 2] <- (1 - v2 * c1)/(v1 * k)
            G.inv[, 3] <- (alpha/v1) * (((v2 * (log(v1) * c1 - 
                c2))/k) - (1 - v2 * c1)/k^2)
            G <- solve(G.inv)
            VC <- (1/n) * G %*% V %*% t(G)
            index <- match(ci.parameter, c("location", "scale", 
                "shape"))
            sd.param <- sqrt(VC[index, index])
        }
    }, mle = {
        neg.ll <- function(theta, x.weird, n.weird) {
            location <- theta[1]
            scale <- theta[2]
            shape <- theta[3]
            con <- 1 - (shape * (x.weird - location))/scale
            if (any(con < 0)) ret.val <- .Machine$double.xmax else {
                y <- -log(con)/shape
                ret.val <- n.weird * log(scale) + (1 - shape) * 
                  sum(y) + sum(exp(-y))
            }
            ret.val
        }
        grad <- function(theta, x.weird, n.weird) {
            location <- theta[1]
            scale <- theta[2]
            shape <- theta[3]
            y <- -log(1 - (shape * (x.weird - location))/scale)/shape
            vec1 <- exp(-y)
            vec2 <- exp(shape * y)
            P <- n.weird - sum(vec1)
            Q <- sum(vec1 * vec2) - (1 - shape) * sum(vec2)
            R <- n.weird - sum(y) + sum(y * vec1)
            g.loc <- Q/scale
            g.scale <- (P + Q)/(scale * shape)
            g.shape <- (R - (P + Q)/shape)/shape
            c(g.loc, g.scale, g.shape)
        }
        if (abs(shape.pwme) < 0.001) shape.pwme <- sign(shape.pwme) * 
            0.001
        dist.params <- nlminb(start = c(location.pwme, scale.pwme, 
            shape.pwme), objective = neg.ll, lower = c(-Inf, 
            .Machine$double.eps, -Inf), upper = c(Inf, Inf, 1), 
            gradient = grad, x.weird = x, n.weird = n)$par
        names(dist.params) <- c("location", "scale", "shape")
        if (ci) {
            if (dist.params["shape"] >= 0.5) warning(paste("Estimated value of the shape parameter", 
                "is greater than or equal to 0.5.  The necessary", 
                "regularity conditions may not hold for this", 
                "confidence interval to be valid."))
            if (information == "expected") {
                scale <- dist.params["scale"]
                shape <- dist.params["shape"]
                c1 <- 1 - shape
                c2 <- gamma(2 - shape)
                s2 <- scale^2
                sh2 <- shape^2
                p <- c1^2 * gamma(1 - 2 * shape)
                q <- c2 * (digamma(c1) - c1/shape)
                M11 <- (n * p)/s2
                M21 <- (n * (p - c2))/(s2 * shape)
                M22 <- (n * (1 - 2 * c2 + p))/(s2 * sh2)
                M31 <- (-n * (q + p/shape))/(scale * shape)
                M32 <- (n * (1 - .Eulers.constant - (1 - c2)/shape - 
                  q - p/shape))/(scale * sh2)
                M33 <- (n * (pi^2/6 + (1 - .Eulers.constant - 
                  1/shape)^2 + (2 * q)/shape + p/sh2))/sh2
                M <- matrix(c(M11, M21, M31, M21, M22, M32, M31, 
                  M32, M33), 3, 3)
                M.inv <- solve(M)
                index <- match(ci.parameter, c("location", "scale", 
                  "shape"))
                sd.param <- sqrt(M.inv[index, index])
            } else {
                grad.hess <- function(theta, x.weird, n.weird) {
                  location <- theta[1]
                  scale <- theta[2]
                  shape <- theta[3]
                  y <- -log(1 - (shape * (x.weird - location))/scale)/shape
                  vec1 <- exp(-y)
                  vec2 <- exp(shape * y)
                  con1 <- sum(vec1)
                  con2 <- sum(vec2)
                  con3 <- sum(vec1 * vec2)
                  con4 <- sum(y * vec1)
                  con5 <- sum(y)
                  con6 <- sum(y * vec1 * vec2)
                  sts <- scale * shape
                  P <- n.weird - con1
                  Q <- con3 - (1 - shape) * con2
                  R <- n.weird - con5 + con4
                  g.location <- Q/scale
                  g.scale <- (P + Q)/sts
                  g.shape <- (R - (P + Q)/shape)/shape
                  pP.loc <- -con3/scale
                  pP.scale <- con1/sts + pP.loc/shape
                  pQ.loc <- (1 - shape)/scale * (sum(vec1 * vec2^2) + 
                    shape * sum(vec2^2))
                  pQ.scale <- pQ.loc/shape - (1 - shape)/sts * 
                    (con3 + shape * con2)
                  pR.loc <- (con2 - con3 + con6)/scale
                  pR.scale <- -(n.weird - con2 + con4 - con1 + 
                    con3 - con6)/sts
                  pR.shape <- (con5 - con4 + sum(y^2 * vec1) - 
                    scale * pR.scale)/shape
                  h11 <- pQ.loc/scale
                  h21 <- (pP.loc + pQ.loc)/sts
                  h22 <- -g.scale/scale + (pP.scale + pQ.scale)/sts
                  h31 <- (pR.loc + scale * h21)/shape
                  h32 <- (pR.scale - g.scale + scale * h22)/shape
                  h33 <- (pR.shape - g.shape + scale * h32)/shape
                  list(gradient = c(g.location, g.scale, g.shape), 
                    hessian = c(h11, h21, h22, h31, h32, h33))
                }
                V.vec <- grad.hess(dist.params, x, n)$hessian
                V.inv <- solve(matrix(V.vec[c(1, 2, 4, 2, 3, 
                  5, 4, 5, 6)], 3, 3))
                index <- match(ci.parameter, c("location", "scale", 
                  "shape"))
                sd.param <- sqrt(V.inv[index, index])
            }
        }
    }, tsoe = {
        if (!is.vector(plot.pos.cons, mode = "numeric") || length(plot.pos.cons) != 
            2) stop("'plot.pos.cons' must be a numeric vector of length 2")
        if (any(is.na(match(c("a", "b"), names(plot.pos.cons))))) names(plot.pos.cons) <- c("a", 
            "b")
        p <- ((1:n) - plot.pos.cons["a"])/(n + plot.pos.cons["b"])
        sort.x <- sort(x)
        dist.params.mat <- matrix(as.numeric(NA), nrow = n - 
            2, ncol = 3)
        for (j in 2:(n - 1)) dist.params.mat[j - 1, ] <- egevd.tsoe.init(x = sort.x[c(1, 
            j, n)], p = p[c(1, j, n)])
        switch(tsoe.method, med = {
            dist.params <- apply(dist.params.mat, 2, median)
            ret.list.method <- paste(method, "based on Median")
        }, lms = {
            location.tsoe <- MASS::lmsreg(rep(1, n - 2), dist.params.mat[, 
                1], intercept = FALSE)$coef
            scale.tsoe <- MASS::lmsreg(rep(1, n - 2), dist.params.mat[, 
                2], intercept = FALSE)$coef
            shape.tsoe <- MASS::lmsreg(rep(1, n - 2), dist.params.mat[, 
                3], intercept = FALSE)$coef
            dist.params <- c(location.tsoe, scale.tsoe, shape.tsoe)
            ret.list.method <- paste(method, "based on Least Median Squares")
        }, lts = {
            location.tsoe <- MASS::ltsreg(rep(1, n - 2), dist.params.mat[, 
                1], intercept = FALSE)$coef
            scale.tsoe <- MASS::ltsreg(rep(1, n - 2), dist.params.mat[, 
                2], intercept = FALSE)$coef
            shape.tsoe <- MASS::ltsreg(rep(1, n - 2), dist.params.mat[, 
                3], intercept = FALSE)$coef
            dist.params <- c(location.tsoe, scale.tsoe, shape.tsoe)
            ret.list.method <- paste(method, "based on Least Trimmed Squares")
        })
        names(dist.params) <- c("location", "scale", "shape")
        if (ci) {
            stop("CI's not yet available when method='tsoe'")
        }
    })
    ret.list <- list(distribution = "Generalized Extreme Value", 
        sample.size = n, parameters = dist.params, n.param.est = 3, 
        method = ret.list.method, data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method)
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        ci.obj <- switch(ci.method, normal.approx = ci.normal.approx(dist.params[ci.parameter], 
            sd.param, n = n, df = n - 2, ci.type = ci.type, alpha = 1 - 
                conf.level))
        ci.obj$parameter <- ci.parameter
        if (method == "mle") 
            ci.obj$method <- paste(ci.obj$method, " based on\n", 
                space(33), information, " information", sep = "")
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
