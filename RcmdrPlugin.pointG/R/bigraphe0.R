
bigraphe0<-function (X) 
{
    p <- ncol(X)
    noms <- colnames(X)
    p.value <- matrix(0, nrow = p, ncol = p)
    colnames(p.value) <- noms
    rownames(p.value) <- noms
    for (j in (1:(p - 1))) {
        for (i in ((j + 1):p)) {
            if (is.factor(X[, i]) & is.factor(X[, j])) {
                ch <- chisq.test(table(X[, i], X[, j]))
                p.value[i, j] <- ch$p.value
                p.value[j, i] <- p.value[i, j]
            }
            if (is.numeric(X[, i]) & is.numeric(X[, j])) {
                fff <- summary(lm(X[, j] ~ X[, i]))$fstatistic
                p.value[i, j] <- pf(q = fff[1], df1 = fff[2], 
                  df2 = fff[3], lower.tail = FALSE)
                p.value[j, i] <- p.value[i, j]
            }
            if (is.numeric(X[, i]) & is.factor(X[, j])) {
                fff <- summary(lm(X[, i] ~ X[, j]))$fstatistic
                p.value[i, j] <- pf(q = fff[1], df1 = fff[2], 
                  df2 = fff[3], lower.tail = FALSE)
                p.value[j, i] <- p.value[i, j]
            }
            if (is.numeric(X[, j]) & is.factor(X[, i])) {
                fff <- summary(lm(X[, j] ~ X[, i]))$fstatistic
                p.value[i, j] <- pf(q = fff[1], df1 = fff[2], 
                  df2 = fff[3], lower.tail = FALSE)
                p.value[j, i] <- p.value[i, j]
            }
        }
    }
    qgraph(1 - p.value, minimum = 0.95, vsize = 10,gray=TRUE,color=rev(brewer.pal(3,name="PuRd"))[1],GLratio=3,labels=noms)
    return(p.value)
}
