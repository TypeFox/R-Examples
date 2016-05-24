QlossRigid <-
function(y1, x0, p1, p0, ...) {

    DX <- y1 - x0

    # sigma2.e <- var(c(DX), na.rm = TRUE)
    # if(sigma2.e == 0) sigma2.e <- 1e-8

    sDX2 <- sum(colSums(DX^2, na.rm = TRUE), na.rm = TRUE)
    N <- sum(colSums(!is.na(DX)))

    # res <- -sDX2/(2 * sigma2.e) - (N / 2) * log(sigma2.e)

    res <- sDX2 / N

    return(res)

}
QcorrRigid <-
function(y1, x0, p1, p0, ...) {

    return(1 - abs(cor(c(y1), c(x0))))

}
