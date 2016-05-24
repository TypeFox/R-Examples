Yule01<-
function (X,Y, levX,nameY) 
{
    CT <- table(X, Y)
    sumX <- apply(CT, 1, sum)
    sumY <- apply(CT, 2, sum)
    sumCT <- sum(CT)
    levY <- levels(Y)
    J <- length(levY)
    result <- matrix(0, nrow = J, ncol = 4)
    for (j in 1:J) {
        a <- CT[levX, levY[j]]
        b <- sumX[levX] - a
        c <- sumY[levY[j]] - a
        d <- sumCT - a - b - c
        result[j, 1] <- a
        Q <- (a * d - b * c)/(a * d + b * c)
        result[j, 2] <- Q
        result[j, 3] <- (1 - Q^2) * sqrt(1/a + 1/b + 1/c + 1/d)/2
        result[j, 4] <- chisq.test(matrix(c(a, c, b, d), ncol = 2))$p.value
    }
    colnames(result) <- c("n", "Q", "se(Q)", "p")
    rownames(result) <- paste(nameY,levY, sep = "_")
    result
}
