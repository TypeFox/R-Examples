`logitreg` <- function(object, groups, k = 1, ...)
    UseMethod("logitreg")

`logitreg.default` <- function(object, groups, k = 1, biasReduced = FALSE,
                               ...) {
    if(!is.factor(groups))
        groups <- factor(groups)
    lev <- levels(groups)
    ## bias reduced fitting via brglm?
    if(biasReduced) {
        FIT <- brglm
    } else {
        FIT <- glm
    }
    within <- without <- vector(mode = "list", length = length(lev))
    names(within) <- names(without) <- lev
    models <- vector(mode = "list", length = length(lev) + 1)
    names(models) <- c(lev, "Combined")
    k <- seq_len(k) + 1
    for(l in lev) {
        inds <- groups == l
        IN <- as.numeric(apply(object[inds, inds], 2,
                               function(x, k) {x[order(x)[k]]}, k = k))
        OUT <- as.numeric(apply(object[inds, !inds], 2,
                                function(x, k) {x[order(x)[k]]}, k = k))
        analogs <- rep(c(TRUE, FALSE), times = c(length(IN), length(OUT)))
        Dij <- c(IN, OUT)
        models[[l]] <- FIT(analogs ~ Dij, data = data.frame(analogs, Dij),
                           family = binomial(link = "logit"), ...)
        models[[l]]$Dij <- Dij
        within[[l]] <- IN
        without[[l]] <- OUT
    }
    IN <- do.call(c, within)
    OUT <- do.call(c, without)
    analogs <- rep(c(TRUE, FALSE), times = c(length(IN), length(OUT)))
    Dij <- c(IN, OUT)
    models[["Combined"]] <- FIT(analogs ~ Dij,
                                data = data.frame(analogs, Dij),
                                family = binomial(link = "logit"), ...)
    models[["Combined"]]$Dij <- Dij
    ##class(models) <- "logitreg"
    out <- list(models = models, groups = groups, method = NULL)
    class(out) <- "logitreg"
    if(!is.null(attr(object, "method"))){
        out$method <- attr(object, "method")
        ## attr(models, "method") <- attr(object, "method")
    }
    out
}

print.logitreg <- function(x, ...) {
    nams <- names(x$models)
    N <- length(x$models)
    cat("\n")
    writeLines(strwrap("Object of class: \"logitreg\""))
    writeLines(strwrap(paste("Number of models:", N)))
    cat("\n")
    writeLines(strwrap("For groups:"))
    print(nams, ...)
    cat("\n")
    invisible(x)
}
