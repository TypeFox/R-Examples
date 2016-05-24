multipleComp<-
function (fit, conf.level = 0.95) 
{
    if (nrow(anova(fit)) != 2) 
        stop("This is not a 1-way ANOVA fit")
    y <- fit$model[, 1]
    f <- fit$model[, 2]
    f <- factor(f)
    fnames <- levels(f)
    k <- length(unique(f))
    df <- length(y) - k
    nr <- k * (k - 1)/2
    contrast.matrix <- matrix(0, nrow = nr, ncol = k)
    row <- 1
    names <- NULL
    for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
            contrast.matrix[row, i] <- 1
            contrast.matrix[row, j] <- -1
            names <- c(names, paste(fnames[i], " - ", fnames[j]))
            row <- row + 1
        }
    }
    row.names(contrast.matrix) <- names
    contrast.matrix <- as.matrix(contrast.matrix)
    estimateContrasts(contrast.matrix, fit, row = TRUE, alpha=1-conf.level)[,1:4]

}

