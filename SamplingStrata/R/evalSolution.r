evalSolution <- function (frame, outstrata, nsampl = 100, writeFiles = FALSE) 
{
    numY <- length(grep("Y", toupper(colnames(frame))))
    numdom <- length(levels(as.factor(frame$DOMAINVALUE)))
    param <- array(0, c(numY, numdom))
    for (i in 1:numY) {
        stmt <- paste("param[i,] <- aggregate(frame$Y", i, ",by=list(frame$DOMAINVALUE),FUN=mean)[,2]", 
            sep = "")
        eval(parse(text = stmt))
    }
    estim <- array(0, c(numdom, nsampl, numY))
    differ <- array(0, c(numdom, nsampl, numY))
    for (j in (1:nsampl)) {
        samp <- selectSample(frame, outstrata, verbatim = FALSE)
        for (k in 1:numY) {
            stmt <- paste("estim[,j,k] <- aggregate(samp$Y", 
                k, "*samp$WEIGHT,by=list(samp$DOMAINVALUE),FUN=sum)[,2] / aggregate(samp$WEIGHT,by=list(samp$DOMAINVALUE),FUN=sum)[,2]", 
                sep = "")
            eval(parse(text = stmt))
        }
    }
    for (j in (1:nsampl)) {
        for (i in (1:numdom)) {
            for (k in 1:numY) {
                differ[i, j, k] <- estim[i, j, k] - param[k, 
                  i]
            }
        }
    }
    cv <- array(0, c(numdom, numY))
    for (k in 1:numY) {
        for (i in (1:numdom)) {
            cv[i, k] <- sd(estim[i, , k])/mean(estim[i, , k])
        }
    }
    cv <- as.data.frame(cv)
    for (i in 1:numY) {
        stmt <- paste("colnames(cv)[", i, "] <- c('CV", i, "')", 
            sep = "")
        eval(parse(text = stmt))
    }
    cv$dom <- paste("DOM", c(1:numdom), sep = "")
    if (writeFiles == TRUE) 
        write.table(cv, "expected_cv.csv", sep = ",", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
    cv1 <- NULL
    cv1$domainvalue <- rep((1:numdom), numY)
    cv1$cv <- rep(0, (numY * numdom))
    cv1$val <- rep(0, (numY * numdom))
    cv1 <- as.data.frame(cv1)
    for (k in (1:numY)) {
        for (j in (1:numdom)) {
            cv1$cv[(k - 1) * numdom + j] <- k
            stmt <- paste("cv1$val[(k-1)*numdom + j] <- cv$CV", 
                k, "[j]", sep = "")
            eval(parse(text = stmt))
        }
    }
    if (writeFiles == TRUE) 
        write.table(cv1, "expected_cv1.csv", sep = ",", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
    diff <- NULL
    diff$dom <- rep(0, numdom * nsampl)
    diff$samp <- rep(0, numdom * nsampl)
    for (i in (1:numY)) {
        stmt <- paste("diff$diff", i, " <- rep(0,numdom*nsampl)", 
            sep = "")
        eval(parse(text = stmt))
    }
    diff <- as.data.frame(diff)
    for (i in (1:numdom)) {
        for (j in (1:nsampl)) {
            diff$dom[(i - 1) * nsampl + j] <- i
            diff$samp[(i - 1) * nsampl + j] <- j
            for (k in (1:numY)) {
                stmt <- paste("diff$diff", k, "[(i-1)*nsampl+j] <- differ[i,j,k]", 
                  sep = "")
                eval(parse(text = stmt))
            }
        }
    }
    if (writeFiles == TRUE) 
        write.table(diff, "differences.csv", sep = ",", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
    if (numdom > 1) {
        if (writeFiles == TRUE) 
            pdf("cv.pdf", width = 7, height = 5)
        boxplot(val ~ cv, data = cv1, col = "orange", main = "Distribution of CV's in the domains", 
            xlab = "Variables Y", ylab = "Value of CV")
        if (writeFiles == TRUE) 
            dev.off()
    }
    if (writeFiles == TRUE) 
        pdf("differences.pdf", width = 14, height = 10)
    k <- ceiling(numY/4)
    for (j in 1:k) {
        split.screen(c(2, 2))
        for (i in 1:4) {
            if (i + 4 * (j - 1) <= numY) {
                stmt <- paste("screen(", i, ")", sep = "")
                eval(parse(text = stmt))
                stmt <- paste("boxplot(diff", i, "~dom,data=diff,ylab='Differences',xlab='Domain',col = 'orange')", 
                  sep = "")
                eval(parse(text = stmt))
                stmt <- paste("mtext(expression(Y", i, "), side=3, adj=0, cex=1.0, line=1)", 
                  sep = "")
                eval(parse(text = stmt))
            }
        }
        close.screen(all.screens = TRUE)
    }
    if (writeFiles == TRUE) 
        dev.off()
    results <- list(coeff_var = cv1, bias = diff)
    return(results)
}
