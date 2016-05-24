# Last modified Feb 16, 2008

`numSummaryIPSUR` <-
function (data, statistics = c("mean", "sd", "skewness", "kurtosis", 
    "quantiles"), quantiles = c(0, 0.25, 0.5, 0.75, 1), groups) 
{
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("abind package missing")
  } else if (!requireNamespace("e1071", quietly = TRUE)){
        stop("e1071 package missing")
  } else {
    data <- as.data.frame(data)
    if (!missing(groups)) 
        groups <- as.factor(groups)
    variables <- names(data)
    statistics <- match.arg(statistics, c("mean", "sd", "skewness", 
        "kurtosis", "quantiles"), several.ok = TRUE)
    ngroups <- if (missing(groups)) 
        1
    else length(grps <- levels(groups))
    quantiles <- if ("quantiles" %in% statistics) 
        quantiles
    else NULL
    quants <- if (length(quantiles) > 1) 
        paste(100 * quantiles, "%", sep = "")
    else NULL
    nquants <- length(quants)
    stats <- c(c("mean", "sd", "skewness", "kurtosis")[c("mean", 
        "sd", "skewness", "kurtosis") %in% statistics], quants)
    nstats <- length(stats)
    nvars <- length(variables)
    result <- list()
    if ((ngroups == 1) && (nvars == 1) && (length(statistics) == 
        1)) {
        if (statistics == "quantiles") 
            table <- quantile(data[, variables], probs = quantiles, 
                na.rm = TRUE)
        else {
            table <- do.call(statistics, list(x = data[, variables], 
                na.rm = TRUE))
            names(table) <- statistics
        }
        NAs <- sum(is.na(data[, variables]))
        n <- nrow(data) - NAs
        result$type <- 1
    }
    else if ((ngroups > 1) && (nvars == 1) && (length(statistics) == 
        1)) {
        if (statistics == "quantiles") {
            table <- matrix(unlist(tapply(data[, variables], 
                groups, quantile, probs = quantiles, na.rm = TRUE)), 
                ngroups, nquants, byrow = TRUE)
            rownames(table) <- grps
            colnames(table) <- quants
        }
        else table <- tapply(data[, variables], groups, statistics, 
            na.rm = TRUE)
        NAs <- tapply(data[, variables], groups, function(x) sum(is.na(x)))
        n <- table(groups) - NAs
        result$type <- 2
    }
    else if ((ngroups == 1)) {
        table <- matrix(0, nvars, nstats)
        rownames(table) <- if (length(variables) > 1) 
            variables
        else ""
        colnames(table) <- stats
        if ("mean" %in% stats) 
            table[, "mean"] <- mean(data[, variables], na.rm = TRUE)
        if ("sd" %in% stats) 
            table[, "sd"] <- sd(data[, variables], na.rm = TRUE)
        if ("skewness" %in% stats) 
            table[, "skewness"] <- apply(as.matrix(data[, variables]), 
                MARGIN = 2, e1071::skewness, na.rm = TRUE)
        if ("kurtosis" %in% stats) 
            table[, "kurtosis"] <- apply(as.matrix(data[, variables]), 
                MARGIN = 2, e1071::kurtosis, na.rm = TRUE)
        if ("quantiles" %in% statistics) {
            table[, quants] <- t(apply(data[, variables, drop = FALSE], 
                2, quantile, probs = quantiles, na.rm = TRUE))
        }
        NAs <- colSums(is.na(data[, variables, drop = FALSE]))
        n <- nrow(data) - NAs
        result$type <- 3
    }
    else {
        table <- array(0, c(ngroups, nstats, nvars), dimnames = list(Group = grps, 
            Statistic = stats, Variable = variables))
        NAs <- matrix(0, nvars, ngroups)
        rownames(NAs) <- variables
        colnames(NAs) <- grps
        for (variable in variables) {
            if ("mean" %in% stats) 
                table[, "mean", variable] <- tapply(data[, variable], 
                  groups, mean, na.rm = TRUE)
            if ("sd" %in% stats) 
                table[, "sd", variable] <- tapply(data[, variable], 
                  groups, sd, na.rm = TRUE)
            if ("skewness" %in% stats) 
                table[, "skewness", variable] <- tapply(data[, 
                  variable], groups, e1071::skewness, na.rm = TRUE)
            if ("kurtosis" %in% stats) 
                table[, "kurtosis", variable] <- tapply(data[, 
                  variable], groups, e1071::kurtosis, na.rm = TRUE)
            if ("quantiles" %in% statistics) {
                res <- matrix(unlist(tapply(data[, variable], 
                  groups, quantile, probs = quantiles, na.rm = TRUE)), 
                  ngroups, nquants, byrow = TRUE)
                table[, quants, variable] <- res
            }
            NAs[variable, ] <- tapply(data[, variable], groups, 
                function(x) sum(is.na(x)))
        }
        if (nstats == 1) 
            table <- table[, 1, ]
        if (nvars == 1) 
            table <- table[, , 1]
        n <- table(groups)
        n <- matrix(n, nrow = nrow(NAs), ncol = ncol(NAs), byrow = TRUE)
        n <- n - NAs
        result$type <- 4
    }
    result$table <- table
    result$statistics <- statistics
    result$n <- n
    if (any(NAs > 0)) 
        result$NAs <- NAs
    class(result) <- "numSummaryIPSUR"
    result

  }
}

`print.numSummaryIPSUR` <-
function (x, ...) 
{
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("abind package missing")
  } else {
    NAs <- x$NAs
    table <- x$table
    n <- x$n
    statistics <- x$statistics
    switch(x$type, "1" = {
        if (!is.null(NAs)) {
            table <- c(table, n, NAs)
            names(table)[length(table) - 1:0] <- c("n", "NA")
        }
        print(table)
    }, "2" = {
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
            rownames(table)[c(1, nrow(table))] <- c(statistics, 
                "n")
            if (!is.null(NAs)) {
                table <- rbind(table, NAs)
                rownames(table)[nrow(table)] <- "NA"
            }
            table <- t(table)
        }
        print(table)
    }, "3" = {
        table <- cbind(table, n)
        colnames(table)[ncol(table)] <- "n"
        if (!is.null(NAs)) {
            table <- cbind(table, NAs)
            colnames(table)[ncol(table)] <- "NA"
        }
        print(table)
    }, "4" = {
        if (length(dim(table)) == 2) {
            table <- cbind(table, t(n))
            colnames(table)[ncol(table)] <- "n"
            if (!is.null(NAs)) {
                table <- cbind(table, t(NAs))
                colnames(table)[ncol(table)] <- "NA"
            }
            print(table)
        }
        else {
            table <- abind::abind(table, t(n), along = 2)
            dimnames(table)[[2]][dim(table)[2]] <- "n"
            if (!is.null(NAs)) {
                table <- abind::abind(table, t(NAs), along = 2)
                dimnames(table)[[2]][dim(table)[2]] <- "NA"
            }
            nms <- dimnames(table)[[3]]
            for (name in nms) {
                cat("\nVariable:", name, "\n")
                print(table[, , name])
            }
        }
    })
    invisible(x)
  }
}
