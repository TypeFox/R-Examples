# various numeric summary statistics

# last modified 2014-09-04 by J. Fox

numSummary <- function(data, 
    statistics=c("mean", "sd", "se(mean)", "IQR", "quantiles", "cv", "skewness", "kurtosis"),
    type=c("2", "1", "3"),
    quantiles=c(0, .25, .5, .75, 1), groups){
    sd <- function(x, type, ...){
        apply(as.matrix(x), 2, stats::sd, na.rm=TRUE)
    }
    IQR <- function(x, type, ...){
        apply(as.matrix(x), 2, stats::IQR, na.rm=TRUE)
    }
    std.err.mean <- function(x, ...){
      x <- as.matrix(x)
      sd <- sd(x)
      n <- colSums(!is.na(x))
      sd/sqrt(n)
    }
    cv <- function(x, ...){
        x <- as.matrix(x)
        mean <- colMeans(x, na.rm=TRUE)
        sd <- sd(x)
        if (any(x <= 0, na.rm=TRUE)) warning("not all values are positive")
        cv <- sd/mean
        cv[mean <= 0] <- NA
        cv
    }
    skewness <- function(x, type, ...){
        if (is.vector(x)) return(e1071::skewness(x, type=type, na.rm=TRUE))
        apply(x, 2, skewness, type=type)
    }
    kurtosis <- function(x, type, ...){
        if (is.vector(x)) return(e1071::kurtosis(x, type=type, na.rm=TRUE))
        apply(x, 2, kurtosis, type=type)
    }
    data <- as.data.frame(data)
    if (!missing(groups)) {
        groups <- as.factor(groups)
        counts <- table(groups)
        if (any(counts == 0)){
            levels <- levels(groups)
            warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
            groups <- factor(groups, levels=levels[counts != 0])
        }
    }
    variables <- names(data)
    if (missing(statistics)) statistics <- c("mean", "sd", "quantiles", "IQR")
    statistics <- match.arg(statistics, c("mean", "sd", "se(mean)", "IQR", "quantiles", "cv", "skewness", "kurtosis"),
        several.ok=TRUE)
    type <- match.arg(type)
    type <- as.numeric(type)
    ngroups <- if(missing(groups)) 1 else length(grps <- levels(groups))
    quantiles <- if ("quantiles" %in% statistics) quantiles else NULL
    quants <- if (length(quantiles) > 1) paste(100*quantiles, "%", sep="") else NULL
    #    quants <- paste(100*quantiles, "%", sep="")
    nquants <- length(quants)
    stats <- c(c("mean", "sd", "se(mean)", "IQR", "cv", "skewness", "kurtosis")[c("mean", "sd", "se(mean)", "IQR", "cv", "skewness", "kurtosis") %in% statistics], quants)
    nstats <- length(stats)
    nvars <- length(variables)
    result <- list()
    if ((ngroups == 1) && (nvars == 1) && (length(statistics) == 1)){
        if (statistics == "quantiles")
            table <- quantile(data[,variables], probs=quantiles, na.rm=TRUE)
        else {
            stats <- statistics
            stats[stats == "se(mean)"] <- "std.err.mean"
            table <- do.call(stats, list(x=data[,variables], na.rm=TRUE, type=type))
            names(table) <- statistics
        }
        NAs <- sum(is.na(data[,variables]))
        n <- nrow(data) - NAs
        result$type <- 1
    }
    else if ((ngroups > 1)  && (nvars == 1) && (length(statistics) == 1)){
        if (statistics == "quantiles"){
            table <- matrix(unlist(tapply(data[, variables], groups,
                quantile, probs=quantiles, na.rm=TRUE)), ngroups, nquants,
                byrow=TRUE)
            rownames(table) <- grps
            colnames(table) <- quants
        }
        else table <- tapply(data[,variables], groups, statistics,
            na.rm=TRUE, type=type)
        NAs <- tapply(data[, variables], groups, function(x)
            sum(is.na(x)))
        n <- table(groups) - NAs
        result$type <- 2
    }
    else if ((ngroups == 1) ){
        X <- as.matrix(data[, variables])
        table <- matrix(0, nvars, nstats)
        rownames(table) <- if (length(variables) > 1) variables else ""
        colnames(table) <- stats
        if ("mean" %in% stats) table[,"mean"] <- colMeans(X, na.rm=TRUE)
        if ("sd" %in% stats) table[,"sd"] <- sd(X)
        if ("se(mean)" %in% stats) table[, "se(mean)"] <- std.err.mean(X)
        if ("IQR" %in% stats) table[, "IQR"] <- IQR(X)
        if ("cv" %in% stats) table[,"cv"] <- cv(X)
        if ("skewness" %in% statistics) table[, "skewness"] <- skewness(X, type=type)
        if ("kurtosis" %in% statistics) table[, "kurtosis"] <- kurtosis(X, type=type)
        if ("quantiles" %in% statistics){
            table[,quants] <- t(apply(data[, variables, drop=FALSE], 2, quantile,
                probs=quantiles, na.rm=TRUE))
        }
        NAs <- colSums(is.na(data[, variables, drop=FALSE]))
        n <- nrow(data) - NAs
        result$type <- 3
    }
    else {
        table <- array(0, c(ngroups, nstats, nvars),
            dimnames=list(Group=grps, Statistic=stats, Variable=variables))
        NAs <- matrix(0, nvars, ngroups)
        rownames(NAs) <- variables
        colnames(NAs) <- grps
        for (variable in variables){
            if ("mean" %in% stats)
                table[, "mean", variable] <- tapply(data[, variable],
                    groups, mean, na.rm=TRUE)
            if ("sd" %in% stats)
                table[, "sd", variable] <- tapply(data[, variable],
                    groups, sd, na.rm=TRUE)
            if ("se(mean)" %in% stats)
              table[, "se(mean)", variable] <- tapply(data[, variable],
                    groups, std.err.mean, na.rm=TRUE)
            if ("IQR" %in% stats)
                table[, "IQR", variable] <- tapply(data[, variable],
                    groups, IQR, na.rm=TRUE)
            if ("cv" %in% stats)
                table[, "cv", variable] <- tapply(data[, variable],
                    groups, cv)
            if ("skewness" %in% stats)
                table[, "skewness", variable] <- tapply(data[, variable],
                    groups, skewness, type=type)
            if ("kurtosis" %in% stats)
                table[, "kurtosis", variable] <- tapply(data[, variable],
                    groups, kurtosis, type=type)
            if ("quantiles" %in% statistics) {
                res <- matrix(unlist(tapply(data[, variable], groups,
                    quantile, probs=quantiles, na.rm=TRUE)), ngroups, nquants,
                    byrow=TRUE)
                table[, quants, variable] <- res
            }
            NAs[variable,] <- tapply(data[, variable], groups, function(x)
                sum(is.na(x)))
        }
        if (nstats == 1) table <- table[,1,]
        if (nvars == 1) table <- table[,,1]
        n <- table(groups)
        n <- matrix(n, nrow=nrow(NAs), ncol=ncol(NAs), byrow=TRUE)
        n <- n - NAs
        result$type <- 4
    }
    result$table <- table
    result$statistics <- statistics
    result$n <- n
    if (any(NAs > 0)) result$NAs <- NAs
    class(result) <- "numSummary"
    result
}

print.numSummary <- function(x, ...){
    NAs <- x$NAs
    table <- x$table
    n <- x$n
    statistics <- x$statistics
    switch(x$type,
        "1" = {
            if (!is.null(NAs)) {
                table <- c(table, n, NAs)
                names(table)[length(table) - 1:0] <- c("n", "NA")
            }
            print(table)
        },
        "2" = {
            if (statistics == "quantiles") {
                table <- cbind(table, n)
                colnames(table)[ncol(table)] <- "n"
                if (!is.null(NAs)) {
                    table <- cbind(table, NAs)
                    colnames(table)[ncol(table)] <- "NA"
                }
            }
            else {
                table <- rbind(table, n)
                rownames(table)[c(1, nrow(table))] <- c(statistics, "n")
                if (!is.null(NAs)) {
                    table <- rbind(table, NAs)
                    rownames(table)[nrow(table)] <- "NA"
                }
                table <- t(table)
            }
            print(table)
        },
        "3" = {
            table <- cbind(table, n)
            colnames(table)[ncol(table)] <- "n"
            if (!is.null(NAs)) {
                table <- cbind(table, NAs)
                colnames(table)[ncol(table)] <- "NA"
            }
            print(table)
        },
        "4" = {
            if (length(dim(table)) == 2){
                n <- t(n)
                nms <- colnames(n)
                colnames(n) <- paste(nms, ":n", sep="")
                table <- cbind(table, n)
                if (!is.null(NAs)) {
                    NAs <- t(NAs)
                    nms <- colnames(NAs)
                    colnames(NAs) <- paste(nms, ":NA", sep="")
                    table <- cbind(table, NAs)
                }
                print(table)
            }
            else {
                table <- abind(table, t(n), along=2)
                dimnames(table)[[2]][dim(table)[2]] <- "n"
                if (!is.null(NAs)) {
                    table <- abind(table, t(NAs), along=2)
                    dimnames(table)[[2]][dim(table)[2]] <- "NA"
                }
                nms <- dimnames(table)[[3]]
                for (name in nms){
                    cat("\nVariable:", name, "\n")
                    print(table[,,name])
                }
            }
        }
    )
    invisible(x)
}
