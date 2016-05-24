summary2way<-function (fit, page = "table", digit = 5, conf.level = 0.95, print.out = TRUE, 
    ...) 
{
    if (!inherits(fit, "lm")) 
        stop("Input is not an \"lm\" object")
    alist <- anova(fit)
    if (nrow(alist) < 3 | nrow(alist) > 4) 
        stop("Not a Two-way ANOVA problem!")
    if (nrow(alist) == 3) {
        if (length(fit$contrast) != 2) 
            stop("All explanatory variables should be factors!")
        inter <- FALSE
    }
    if (nrow(alist) == 4) {
        if (alist$Df[4] == 0) 
            stop("Need more than one observation per cell for interaction!")
        if (length(fit$model) != 3) 
            stop("Not a Two-way ANOVA problem!")
        if (length(fit$contrast) != 2) 
            stop("All explanatory variables should be factors!")
        inter <- TRUE
    }
    y <- fit$model[, 1]
    f1 <- factor(fit$model[, 2])
    f2 <- factor(fit$model[, 3])
    f1f2 <- as.factor(crossFactors(f1, f2))
    nlevf1 <- length(unique(f1))
    nlevf2 <- length(unique(f2))
    if (inter) 
        m <- 3
    else m <- 2
    a.df <- c(alist$Df, sum(alist$Df))
    a.ss <- round(c(alist$"Sum Sq", sum(alist$"Sum Sq")), digit)
    a.ms <- round(alist$"Mean Sq", digit)
    fvalue <- round(alist$"F value"[1:m], digit)
    pvalue <- round(alist$"Pr(>F)"[1:m], digit)
    a.table <- cbind(a.df, a.ss, c(paste(a.ms), ""), c(paste(fvalue), 
        "", ""), c(paste(pvalue), "", ""))
    if (inter) 
        dimnames(a.table) <- list(c(row.names(alist), "Total        "), 
            c("Df ", "Sum Squares ", "Mean Square ", "F-statistic ", 
                "p-value   "))
    else dimnames(a.table) <- list(c(row.names(alist), "Total    "), 
        c("Df ", "Sum Squares ", "Mean Square ", "F-statistic ", 
            "p-value   "))
    group1 <- split(y, f1)
    group2 <- split(y, f2)
    n1 <- length(group1)
    n2 <- length(group2)
    if (inter) {
        f <- factor(crossFactors(f1, f2))
        group3 <- split(y, f)
        n3 <- length(group3)
        group <- c(group1, group2, group3)
    }
    else group <- c(group1, group2)
    n <- length(group)
    size <- c(length(y), numeric(n))
    mea <- c(mean(y), numeric(n))
    med <- c(median(y), numeric(n))
    std <- c(sd(y), numeric(n))
    mid <- c(quantile(y, 0.75) - quantile(y, 0.25), numeric(n))
    for (i in 2:(n + 1)) {
        size[i] <- length(group[[i - 1]])
        g <- numeric(size[i])
        g <- group[[i - 1]]
        mea[i] <- mean(g)
        med[i] <- median(g)
        std[i] <- sd(g)
        mid[i] <- quantile(g, 0.75) - quantile(g, 0.25)
    }
    size <- round(size, digit)
    mea <- round(mea, digit)
    med <- round(med, digit)
    std <- round(std, digit)
    mid <- round(mid, digit)
    if (inter) {
        numeric.summary <- cbind(c(size[1], "", size[2:(n1 + 
            1)], "", size[(2 + n1):(n1 + n2 + 1)], "", size[(2 + 
            n1 + n2):(n + 1)]), c(mea[1], "", mea[2:(n1 + 1)], 
            "", mea[(2 + n1):(n1 + n2 + 1)], "", mea[(2 + n1 + 
                n2):(n + 1)]), c(med[1], "", med[2:(n1 + 1)], 
            "", med[(2 + n1):(n1 + n2 + 1)], "", med[(2 + n1 + 
                n2):(n + 1)]), c(std[1], "", std[2:(n1 + 1)], 
            "", std[(2 + n1):(n1 + n2 + 1)], "", std[(2 + n1 + 
                n2):(n + 1)]), c(mid[1], "", mid[2:(n1 + 1)], 
            "", mid[(2 + n1):(n1 + n2 + 1)], "", mid[(2 + n1 + 
                n2):(n + 1)]))
        dimnames(numeric.summary) <- list(c("All Data  ", paste("By ", 
            attributes(fit$terms)$variables[[3]], ":  "), paste(names(group1), 
            "  "), paste("By ", attributes(fit$terms)$variables[[4]], 
            ":  "), paste(names(group2), "  "), "Combinations:  ", 
            paste(names(group3), "  ")), c("Size  ", "Mean  ", 
            "Median  ", "Std Dev  ", "Midspread  "))
    }
    else {
        numeric.summary <- cbind(c(size[1], "", size[2:(n1 + 
            1)], "", size[(2 + n1):(n1 + n2 + 1)]), c(mea[1], 
            "", mea[2:(n1 + 1)], "", mea[(2 + n1):(n1 + n2 + 
                1)]), c(med[1], "", med[2:(n1 + 1)], "", med[(2 + 
            n1):(n1 + n2 + 1)]), c(std[1], "", std[2:(n1 + 1)], 
            "", std[(2 + n1):(n1 + n2 + 1)]), c(mid[1], "", mid[2:(n1 + 
            1)], "", mid[(2 + n1):(n1 + n2 + 1)]))
        dimnames(numeric.summary) <- list(c("All Data  ", paste("By ", 
            attributes(fit$terms)$variables[[3]], ":  "), paste(names(group1), 
            "  "), paste("By ", attributes(fit$terms)$variables[[4]], 
            ":  "), paste(names(group2), "  ")), c("Size  ", 
            "Mean  ", "Median  ", "Std Dev  ", "Midspread  "))
    }
    dc <- dummy.coef(fit)
    if (!inter) 
        dc[[4]] <- rep(0, n1 * n2)
    effmat <- (matrix(dc[[4]], n1, n2) + outer(dc[[2]], rep(1, 
        n2)) + outer(rep(1, n1), dc[[3]])) + dc[[1]]
    mmean <- mean(effmat)
    cellmns <- matrix(NA, length(levels(f1)), length(levels(f2)))
    levf1 <- levels(f1)
    levf2 <- levels(f2)
    for (i in 1:length(levf1)) {
        for (j in 1:length(levf2)) {
            cellmns[i, j] <- mean(y[f1 == levf1[i] & f2 == levf2[j]])
        }
    }
    cellmns <- cbind(cellmns, apply(effmat, 1, mean))
    cellmns <- rbind(cellmns, c(apply(effmat, 2, mean), mmean))
    dimnames(cellmns) <- list(c(rep(" ", length(levels(f1)) + 
        1)), c(rep(" ", length(levels(f2)) + 1)))
    matr <- rbind(rep("", ncol(effmat)))
    matr[1, as.integer((length(levels(f2)) + 1)/2)] <- as.character(attributes(fit$terms)$variables[[4]])
    matr <- c("", "", matr, "")
    matc <- cbind(rep("", nrow(effmat)))
    matc[as.integer((length(levels(f1)) + 1)/2), 1] <- as.character(attributes(fit$terms)$variables[[3]])
    matc <- c(matc, "")
    C <- c(levels(f1), as.character(attributes(fit$terms)$variables[[4]]))
    R <- c("", "", levels(f2), as.character(attributes(fit$terms)$variables[[3]]))
    M1 <- cbind(matc, C, format(round(cellmns, digit), digit = 5))
    M2 <- rbind(matr, R, M1)
    roweff <- apply(effmat, 1, mean) - mean(apply(effmat, 1, 
        mean))
    coleff <- apply(effmat, 2, mean) - mean(apply(effmat, 2, 
        mean))
    interact <- (effmat - outer(apply(effmat, 1, mean), rep(1, 
        n2)) - outer(rep(1, n1), apply(effmat, 2, mean))) + mean(effmat)
    effmat1 <- cbind(interact, roweff)
    effmat2 <- rbind(effmat1, c(coleff, mmean))
    dimnames(effmat2) <- list(c(rep(" ", length(levels(f1)) + 
        1)), c(rep(" ", length(levels(f2)) + 1)))
    matr <- rbind(rep("", ncol(effmat)))
    matr[1, as.integer((length(levels(f2)) + 1)/2)] <- as.character(attributes(fit$terms)$variables[[4]])
    matr <- c("", "", matr, "")
    matc <- cbind(rep("", nrow(effmat)))
    matc[as.integer((length(levels(f1)) + 1)/2), 1] <- as.character(attributes(fit$terms)$variables[[3]])
    matc <- c(matc, "")
    C <- c(levels(f1), paste(as.character(attributes(fit$terms)$variables[[4]]), 
        "effect"))
    R <- c("", "", levels(f2), paste(as.character(attributes(fit$terms)$variables[[3]]), 
        "effect"))
    M3 <- cbind(matc, C, format(round(effmat2, digit), digit = 5))
    M4 <- rbind(matr, R, M3)
    f1.name <- row.names(alist)[1]
    f2.name <- row.names(alist)[2]
    if (page == "table") {
        cat("ANOVA Table:\n")
        print(a.table, quote = FALSE)
    }
    if (page == "means") {
        cat("\n\nCell-means Matrix:\n")
        print(M2, quote = FALSE)
        cat("\nNumeric Summary:\n")
        print(numeric.summary, quote = FALSE)
    }
    if (page == "effects") {
        cat("\n\nTable of Effects:\n")
        print(M4, quote = FALSE)
    }
    if (page == "interaction") {
        cat(paste("\n\nComparisons within ", f1.name, ":\n\n", 
            sep = ""))
        contrast.matrix1 <- names <- NULL
        offset <- 1
        for (levs in 1:nlevf1) {
            temp <- matrix(0, nrow = (nlevf2 * (nlevf2 - 1)/2), 
                ncol = nlevf1 * nlevf2)
            row <- 1
            for (i in offset:(levs * nlevf2 - 1)) {
                for (j in (i + 1):(levs * nlevf2)) {
                  temp[row, i] <- 1
                  temp[row, j] <- -1
                  names <- c(names, paste(levels(f1f2)[i], " - ", 
                    levels(f1f2)[j]))
                  row <- row + 1
                }
            }
            contrast.matrix1 <- as.matrix(rbind(contrast.matrix1, 
                temp))
            offset <- offset + nlevf2
        }
        row.names(contrast.matrix1) <- names
        fit.1way <- lm(y ~ f1f2)
        L <- (nlevf1 * nlevf2/2) * (1 + nlevf1)
        contrasts1 <- estimateContrasts(as.matrix(contrast.matrix1), 
            fit.1way, alpha=1-conf.level,row = TRUE, L)
        print(contrasts1, quote = FALSE)
        cat(paste("\n\nComparisons between ", f1.name, ":\n\n", 
            sep = ""))
        contrast.matrix2 <- names <- NULL
        temp <- matrix(0, nrow = nlevf2, ncol = nlevf1 * nlevf2)
        nrows <- nlevf1 * (nlevf1 - 1)/2
        for (i in seq(1, nlevf1 * nlevf2 - nlevf2, nlevf2)) {
            for (j in seq(i + nlevf2, nlevf1 * nlevf2, nlevf2)) {
                for (row in 1:nlevf2) {
                  temp[row, i + row - 1] <- 1
                  temp[row, j + row - 1] <- -1
                  names <- c(names, paste(levels(f1f2)[i + row - 
                    1], " - ", levels(f1f2)[j + row - 1]))
                }
                contrast.matrix2 <- as.matrix(rbind(contrast.matrix2, 
                  temp))
                temp <- matrix(0, nrow = nlevf2, ncol = nlevf1 * 
                  nlevf2)
            }
        }
        row.names(contrast.matrix2) <- names
        contrasts2 <- estimateContrasts(as.matrix(contrast.matrix2), 
            fit.1way, alpha=1-conf.level,row = TRUE,L)
        print(contrasts2, quote = FALSE)
    }
    if (page == "nointeraction") {
        cat(paste("\n\n", f1.name, " comparisons:\n\n", sep = ""))
        contrast.matrix1 <- matrix(0, nrow = nlevf1 * (nlevf1 - 
            1)/2, ncol = nlevf1)
        row <- 1
        names <- NULL
        for (i in 1:(nlevf1 - 1)) {
            for (j in (i + 1):nlevf1) {
                contrast.matrix1[row, i] <- 1
                contrast.matrix1[row, j] <- -1
                names <- c(names, paste(levels(f1)[i], " - ", 
                  levels(f1)[j]))
                row <- row + 1
            }
        }
        row.names(contrast.matrix1) <- names
        contrast.matrix1 <- as.matrix(contrast.matrix1)
        contrasts1 <- estimateContrasts(contrast.matrix1, fit, alpha=1-conf.level,
            row = TRUE)
        print(contrasts1, quote = FALSE)
        cat(paste("\n\n", f2.name, " comparisons:\n\n", sep = ""))
        contrast.matrix2 <- matrix(0, nrow = nlevf2 * (nlevf2 - 
            1)/2, ncol = nlevf2)
        row <- 1
        names <- NULL
        for (i in 1:(nlevf2 - 1)) {
            for (j in (i + 1):nlevf2) {
                contrast.matrix2[row, i] <- 1
                contrast.matrix2[row, j] <- -1
                names <- c(names, paste(levels(f2)[i], " - ", 
                  levels(f2)[j]))
                row <- row + 1
            }
        }

        row.names(contrast.matrix2) <- names
        contrast.matrix2 <- as.matrix(contrast.matrix2)
        contrasts2 <- estimateContrasts(contrast.matrix2, fit, alpha=1-conf.level,
            row = FALSE)
         
     
        print(contrasts2, quote = FALSE)
    }
    if (!inter) 
        interact <- NULL
    invisible(list(Df = a.df, "Sum of Sq" = a.ss, "Mean Sq" = a.ms, 
        "F value" = alist$"F value"[1:m], "Pr(F)" = alist$"Pr(>F)"[1:m], 
        "Main Effects" = mmean, "Row Effects" = roweff, "Col Effects" = coleff, 
        Interaction = interact))

}

