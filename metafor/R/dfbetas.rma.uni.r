dfbetas.rma.uni <-
function (model, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    tau2.del <- rep(NA_real_, x$k.f)
    dfbetas <- matrix(NA_real_, nrow = x$k.f, ncol = x$p)
    if (x$weighted && !is.null(x$weights)) {
        A <- diag(x$weights, nrow = x$k, ncol = x$k)
        stXAX <- .invcalc(X = x$X, W = A, k = x$k)
    }
    else {
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
    }
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(suppressWarnings(rma(x$yi.f[-i], x$vi.f[-i], 
            weights = x$weights.f[-i], mods = cbind(x$X.f[-i, 
                ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control)), 
            silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        if (any(res$coef.na)) 
            next
        tau2.del[i] <- res$tau2
        dfbeta <- x$b - res$b
        if (x$weighted) {
            if (is.null(x$weights)) {
                vb.del <- .invcalc(X = x$X, W = diag(1/(x$vi + 
                  tau2.del[i]), nrow = x$k, ncol = x$k), k = x$k)
            }
            else {
                vb.del <- tcrossprod(stXAX, x$X) %*% A %*% diag(x$vi + 
                  tau2.del[i], nrow = x$k, ncol = x$k) %*% A %*% 
                  x$X %*% stXAX
            }
        }
        else {
            vb.del <- tcrossprod(stXX, x$X) %*% diag(x$vi + tau2.del[i], 
                nrow = x$k, ncol = x$k) %*% x$X %*% stXX
        }
        dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
    }
    if (na.act == "na.omit") {
        out <- dfbetas[x$not.na, , drop = FALSE]
        rownames(out) <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- dfbetas
        rownames(out) <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    colnames(out) <- rownames(x$b)
    out <- data.frame(out)
    return(out)
}
