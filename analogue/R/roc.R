#######################################################################
## Calculation of the ROC curve points themselves,
## based on slightly modified code of Thomas Lumley
## (used with permission), in his Programmer's Niche article in R News
## (Vol. 4(1) 33--36). Uses the optimisations in the article to
## calculate the ROC curve itself.
#######################################################################

`roc` <- function(object, groups, k = 1, ...) UseMethod("roc")

`roc.default` <- function(object, groups, k = 1,
                          thin = FALSE, max.len = 10000, ...) {
    calcROC <- function(IN, OUT, thin = FALSE, max.len = 10000) {
        n.IN <- length(IN)
        n.OUT <- length(OUT)
        if(thin) {
            ## thin
            ratio <- n.IN/n.OUT
            largest <- which.max(c(n.IN, n.OUT))
            seq.vars <- c(max.len * ratio, max.len)
            if(largest == 1)
                seq.vars <- rev(seq.vars)
            IN <- quantile(IN,
                           probs = seq(0, 1, length.out = seq.vars[1]))
            OUT <- quantile(OUT,
                            probs = seq(0, 1, length.out = 1/seq.vars[2]))
            n.IN <- length(IN)
            n.OUT <- length(OUT)
        }
        g <- rep(c(TRUE, FALSE), times = c(n.IN, n.OUT))
        tab <- table(c(IN, OUT), g)
        TPF <- cumsum(tab[, 2])/sum(tab[, 2])
        FPE <- cumsum(tab[, 1])/sum(tab[, 1])
        roc.values <- TPF - FPE
        ##roc.points <- rev(sort(c(IN, OUT)))
        roc.points <- rev(sort(as.numeric(dimnames(tab)[[1]])))
        optimal <- as.numeric(names(max.roc <- which.max(roc.values)))
        names(FPE) <- names(TPF) <- names(roc.values) <- NULL
        wilcox <- wilcox.test(IN, OUT, conf.int = FALSE)
        AUC <- 1 - (wilcox$statistic / (n.IN * n.OUT))
        AUC2 <- AUC^2
        q1 <- AUC / (2-AUC)
        q2 <- (2 * AUC2) / (1 + AUC)
        se.fit <- AUC * (1 - AUC) + ((n.IN - 1) * (q1 - AUC2)) +
            ((n.OUT - 1) * (q2 - AUC2))
        se.fit <- sqrt(se.fit / (n.IN * n.OUT))
        p.value <- wilcox$p.value
        prior <- c(n.IN, n.OUT) / sum(n.IN, n.OUT)
        retval <- list(TPF = TPF, FPE = FPE, optimal = optimal,
                       AUC = AUC, se.fit = se.fit, n.in = n.IN,
                       n.out = n.OUT, p.value = p.value,
                       roc.points = roc.points, max.roc = max.roc,
                       prior = prior,
                       analogue = list(yes = IN, no = OUT))
        retval
    }
    if(!is.factor(groups))
        groups <- factor(groups)
    lev <- levels(groups)
    n.g <- length(lev)
    within <- without <- vector(mode = "list", length = length(lev))
    names(within) <- names(without) <- lev
    roc <- vector(mode = "list", length = length(lev) + 1)
    names(roc) <- c(lev, "Combined")
    statistics <- data.frame(matrix(NA, nrow = n.g+1, ncol = 6))
    names(statistics) <- c("In","Out","Opt. Dis.","AUC","SE","p-value")
    rownames(statistics) <- c(lev, "Combined")
    k <- seq_len(k) + 1
    for(l in lev) {
        inds <- groups == l
        IN <- as.numeric(apply(object[inds, inds], 2,
                               function(x, k) {x[order(x)[k]]}, k = k))
        OUT <- as.numeric(apply(object[inds, !inds], 2,
                                function(x, k) {x[order(x)[k]]}, k = k))
        ROC <- calcROC(IN, OUT, thin = thin, max.len = max.len)
        within[[l]] <- ROC$analogue$yes
        without[[l]] <- ROC$analogue$no
        statistics[l, ] <- with(ROC, data.frame(n.in, n.out, optimal, AUC,
                                                se.fit, p.value))
        roc[[l]] <- ROC
    }
    IN <- do.call(c, within)
    OUT <- do.call(c, without)
    roc[["Combined"]] <- ROC <- calcROC(IN, OUT, thin = thin,
                                        max.len = max.len)
    statistics["Combined", ] <- with(ROC, data.frame(n.in, n.out, optimal,
                                                     AUC, se.fit, p.value))
    retval <- list(statistics = statistics, roc = roc)
    class(retval) <- "roc"
    if(!is.null(attr(object, "method")))
        attr(retval, "method") <- attr(object, "method")
    return(retval)
}

`roc.mat` <- function(object, groups, k = 1, ...) {
    retval <- roc(object$Dij, groups = groups, k = k, ...)
    attr(retval, "method") <- attr(object, "method")
    return(retval)
}

`roc.analog` <- function(object, groups, k = 1, ...) {
    if(is.null(object$train))
        stop("'object$train' missing. Refit 'object' with argument 'keep.train = TRUE'")
    retval <- roc(object$train, groups = groups, k = k, ...)
    attr(retval, "method") <- object$method
    return(retval)
}
