"print.relimp" <-
function (x, digits = 3, ...)
{
    object <- x
    sets <- object$sets
    labels <- names(sets)
    l1 <- length(sets[[1]])
    l2 <- length(sets[[2]])
    if (l1 > l2) {
        sets[[2]] <- c(sets[[2]], rep("", l1 - l2))
    }
    if (l2 > l1) {
        sets[[1]] <- c(sets[[1]], rep("", l2 - l1))
    }
    sets <- as.data.frame(sets)
    response.cat <- object$response.category
    names(sets) <- c(paste("     Numerator effects (\"", labels[1],
                           "\")", sep = ""),
                     paste("     Denominator effects (\"",
                           labels[2], "\")", sep = ""))
    zz <- capture.output(print(sets, rowlab = rep("", nrow(sets))))
    CI95 <- object$log.ratio +
        1.96 * (object$se.log.ratio) * c(-1, 1)
    cat(paste("\n", "Relative importance summary for model\n",
              "    ", paste(deparse(object$model), collapse = "\n"), "\n",
              if (!is.null(response.cat)) {
                  paste("response category", response.cat, "\n\n")
              },
              if (!is.null(object$dispersion))
              paste("with dispersion set to", object$dispersion,
                    "\n\n")
              else "\n", paste(zz, "", collapse = "\n"),
              "\n\n", "Ratio of effect standard deviations: ",
              round(exp(object$log.ratio), digits),
              "\n", "Log(sd ratio):                 ",
              round(object$log.ratio, digits), "   (se ",
              round(object$se.log.ratio, digits), ")\n\n",
              "Approximate 95% confidence interval for log(sd ratio): (",
              paste(round(CI95, digits), collapse = ","), ")\n",
              "Approximate 95% confidence interval for sd ratio:      (",
              paste(round(exp(CI95), digits), collapse = ","), ")\n",
              sep = "")
        )
}

"relimp" <-
function (object, set1 = NULL, set2 = NULL, label1 = "set1",
    label2 = "set2", subset = TRUE, response.cat = NULL, ...)
{
    if (inherits(object, "multinom")) {
        if (requireNamespace("nnet", quietly = TRUE)) {}
        else
        {stop("the necessary \"nnet\" package is not available")}
    }
    if (inherits(object, "multinom") && is.null(response.cat))
        stop("argument `response.cat' must be specified")
    if (!is.null(response.cat))
        response.cat <- as.character(response.cat) ## numbers get coerced
    if (inherits(object, "multinom") && !(response.cat %in%
                                          rownames(coef(object))))
        stop("argument `response.cat' not valid for this model")
    covmat <- vcov(object)
    if (is.null(set1) || is.null(set2)) {
        coefnames <- {
            if (!inherits(object, "multinom"))
                colnames(covmat)
            else colnames(coef(object))
        }
        if (is.null(coefnames)) coefnames <- names(coef(object))
        sets <- pickFrom(coefnames, nsets = 2, return.indices = FALSE,
                         setlabels = c(label1, label2),
                         title = "Specify a relative importance (\"relimp\") comparison",
                         items.label = "Model coefficients")
        if (!is.null(sets)) {
            set1 <- sets[[1]]
            set2 <- sets[[2]]
            label1 <- names(sets)[1]
            label2 <- names(sets)[2]
        }
        else stop("\neffects for comparison (set1,set2) not specified")
    }
    else {
        if (max(union(set1, set2)) > length(coef(object))) {
            stop("Index out of bounds")
        }
        if (any(is.na(coef(object)[set1]))){
            stop("set1 contains NAs")}
        if (any(is.na(coef(object)[set2]))){
            stop("set2 contains NAs")}
        coefnames <- if (!inherits(object, "multinom")) names(coef(object))
        else colnames(coef(object))
        if (is.null(coefnames)) coefnames <- names(coef(object))
        set1 <- coefnames[set1]
        set2 <- coefnames[set2]}
    ## notation below follows Silber, Rosenbaum and Ross (1995, JASA)
    coefs <- {
        if (!inherits(object, "multinom"))
            coef(object)
        else coef(object)[response.cat, ]
    }
    if (inherits(object, "multinom")) {
        the.coefs <- coef(object)
        indices <- t(matrix(1:prod(dim(the.coefs)),
                            nrow = ncol(the.coefs),
                            ncol = nrow(the.coefs),
                            dimnames = dimnames(t(the.coefs))))
        indices <- indices[response.cat, ]
        covmat <- covmat[indices, indices]
    }
    beta <- coefs[set1, drop = FALSE]
    gamma <- coefs[set2, drop = FALSE]
    if (!is.matrix(object$x))
        modelmatrix <- model.matrix(object)
    else modelmatrix <- object$x
    if (is.numeric(subset) || (is.logical(subset) && (length(subset) ==
                               1 || length(subset) == nrow(modelmatrix)))) {
        X <- modelmatrix[subset, set1, drop = FALSE]
        H <- modelmatrix[subset, set2, drop = FALSE]
    }
    else stop(
     paste("\nspecified subset should be either a vector of numeric indices,",
        "\nor a logical vector with length equal to the number of rows",
        "\nin the model frame for model", paste("\"",
                                                deparse(match.call()$object),
                                                "\"", sep = "")))
    X <- sweep(X, 2, apply(X, 2, mean))
    H <- sweep(H, 2, apply(H, 2, mean))
    indices <- if (inherits(object, "multinom")){
        c(paste(response.cat, set1, sep=":"),
          paste(response.cat, set2, sep=":"))}
    else c(set1, set2)
    Sigma <- covmat[indices, indices]
    pi <- X %*% beta
    phi <- H %*% gamma
    sd.ratio <- sd(as.vector(pi))/sd(as.vector(phi))
    log.ratio <- log(sd.ratio)
    w <- rbind((t(X) %*% pi)/sum(pi * pi), (-t(H) %*% phi)/sum(phi * phi))
    var.log.ratio <- t(w) %*% Sigma %*% w
    se.log.ratio <- sqrt(var.log.ratio)
    Call <- match.call() ## to see if the dispersion parameter was given
    dispersion <- {
        if (pmatch("disp", names(Call), 0) > 0)
            Call$dispersion
        else NULL
    }
    ans <- list(model = object$call, response.category = response.cat,
                dispersion = dispersion, sets = list(set1, set2),
                log.ratio = log.ratio, se.log.ratio = se.log.ratio)
    names(ans$sets) <- c(label1, label2)
    class(ans) <- "relimp"
    return(ans)
}

"relrelimp" <-
function (object, set1 = NULL, set2 = NULL, label1 = "set1",
          label2 = "set2", subset = TRUE,
          response.cat1 = NULL, response.cat2 = NULL)
{
    if (!inherits(object, "multinom")) {
        stop("Object is not of class \"multinom\"")
    }
    if (requireNamespace("nnet", quietly = TRUE)) {}
        else {stop("the \"nnet\" package is not available")}
    if (is.null(response.cat1) || is.null(response.cat2))
        stop("arguments `response.cat1' and `response.cat2' must be specified")
    response.cat1 <- as.character(response.cat1) ## numbers get coerced
    response.cat2 <- as.character(response.cat2)
    the.coefs <- coef(object)
    if (!(response.cat1 %in% rownames(the.coefs)) ||
        !(response.cat2 %in% rownames(the.coefs)))
        stop("`response.cat' argument(s) not valid for this model")
    if (is.null(set1) || is.null(set2)) {
        coefnames <- {
            if (!inherits(object, "multinom"))
                names(coef(object))
            else colnames(the.coefs)
        }
        sets <- pickFrom(coefnames, nsets = 2, return.indices = TRUE,
                         setlabels = c(label1, label2),
                         title = "Specify a relative importance (\"relimp\") comparison",
            items.label = "Model coefficients")
        if (!is.null(sets)) {
            set1 <- sets[[1]]
            set2 <- sets[[2]]
            label1 <- names(sets)[1]
            label2 <- names(sets)[2]
        }
        else stop("\neffects for comparison (set1,set2) not specified")
    }
    if (max(union(set1, set2)) > length(the.coefs)) {
        stop("Index out of bounds")
    }
    ## notation below follows Silber, Rosenbaum and Ross (1995, JASA)
    coefs <- the.coefs[c(response.cat1, response.cat2), ]
    covmat <- vcov(object)
    indices <- t(matrix(1:prod(dim(the.coefs)),
                        nrow = ncol(the.coefs),
                        ncol = nrow(the.coefs),
                        dimnames = dimnames(t(the.coefs))))
    indices <- as.vector(t(indices[c(response.cat1, response.cat2),
                                   ]))
    covmat <- covmat[indices, indices]
    beta1 <- as.vector(coefs[1, set1, drop = FALSE])
    gamma1 <- as.vector(coefs[1, set2, drop = FALSE])
    beta2 <- as.vector(coefs[2, set1, drop = FALSE])
    gamma2 <- as.vector(coefs[2, set2, drop = FALSE])
    names <- colnames(coefs)
    names1 <- names[set1]
    names2 <- names[set2]
    if (!is.matrix(object$x))
        modelmatrix <- model.matrix(object)
    else modelmatrix <- object$x
    if (is.numeric(subset) || (is.logical(subset) &&
                               (length(subset) == 1 ||
                                length(subset) == nrow(modelmatrix)))) {
        X <- modelmatrix[subset, set1, drop = FALSE]
        H <- modelmatrix[subset, set2, drop = FALSE]
    }
    else stop(
              paste("\nspecified subset should be either a vector of numeric indices,",
                    "\nor a logical vector with length equal to the number of rows",
                    "\nin the model frame for model",
                    paste("\"", deparse(match.call()$object),
                          "\"", sep = "")))
    X <- sweep(X, 2, apply(X, 2, mean))
    H <- sweep(H, 2, apply(H, 2, mean))
    set1a <- set1 + ncol(coefs)
    set2a <- set2 + ncol(coefs)
    Sigma <- covmat[c(set1, set2, set1a, set2a),
                    c(set1, set2, set1a, set2a)]
    Sigma1 <- covmat[c(set1, set2), c(set1, set2)]
    Sigma2 <- covmat[c(set1a, set2a), c(set1a, set2a)]
    pi1 <- X %*% beta1
    phi1 <- H %*% gamma1
    pi2 <- X %*% beta2
    phi2 <- H %*% gamma2
    sd.ratio1 <- sd(as.vector(pi1))/sd(as.vector(phi1))
    sd.ratio2 <- sd(as.vector(pi2))/sd(as.vector(phi2))
    log.ratio1 <- log(sd.ratio1)
    log.ratio2 <- log(sd.ratio2)
    log.ratioratio <- log.ratio1 - log.ratio2
    w1 <- rbind((t(X) %*% pi1)/sum(pi1 * pi1),
                (-t(H) %*% phi1)/sum(phi1 * phi1))
    var.log.ratio1 <- t(w1) %*% Sigma1 %*% w1
    se.log.ratio1 <- sqrt(var.log.ratio1)
    w2 <- rbind((t(X) %*% pi2)/sum(pi2 * pi2),
                (-t(H) %*% phi2)/sum(phi2 * phi2))
    var.log.ratio2 <- t(w2) %*% Sigma2 %*% w2
    se.log.ratio2 <- sqrt(var.log.ratio2)
    w <- rbind((t(X) %*% pi1)/sum(pi1 * pi1),
               (-t(H) %*% phi1)/sum(phi1 * phi1),
               (-t(X) %*% pi2)/sum(pi2 * pi2),
               (t(H) %*% phi2)/sum(phi2 * phi2))
    var.log.ratioratio <- t(w) %*% Sigma %*% w
    se.log.ratioratio <- sqrt(var.log.ratioratio)
    Call <- match.call()  ## to see if the dispersion parameter was given
    dispersion <- {
        if (pmatch("disp", names(Call), 0) > 0)
            Call$dispersion
        else NULL
    }
    ans <- list(model = object$call,
                response.category = c(response.cat1, response.cat2),
                dispersion = dispersion,
                sets = list(names1, names2),
                log.ratio = c(log.ratio1, log.ratio2, log.ratioratio),
                se.log.ratio = c(se.log.ratio1,
                                 se.log.ratio2,
                                 se.log.ratioratio)
                )
    names(ans$sets) <- c(label1, label2)
    class(ans) <- "relrelimp"
    return(ans)
}
