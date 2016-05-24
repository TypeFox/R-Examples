variogram <- function (indv, time, Y) 
{
    id <- as.vector(indv[!is.na(Y)])
    time <- as.vector(time[!is.na(Y)])
    y <- as.vector(Y[!is.na(Y)])
    subject <- unique(id)
    m <- length(subject)
    vv <- c()
    vt <- c()
    vtot <- c()
    for (i in 1:m) {
        j1 <- id == subject[i]
        rr1 <- y[j1]
        tt <- time[j1]
        dr <- outer(rr1, rr1, function(x, y) {
            0.5 * (x - y)^2
        })
        vv <- c(vv, dr[upper.tri(dr)])
        dt <- outer(tt, tt, function(x, y) {
            abs(x - y)
        })
        vt <- c(vt, dt[upper.tri(dt)])
        l <- i + 1
        while (l <= m) {
            j2 <- id == subject[l]
            rr2 <- y[j2]
            dtot <- outer(rr1, rr2, function(x, y) {
                0.5 * (x - y)^2
            })
            vtot <- c(vtot, c(dtot))
            l <- l + 1
        }
    }
    svar <- cbind(vt, vv)
    sigma2 <- mean(vtot)
    vrgm <- list(svar = svar, sigma2 = sigma2)
    class(vrgm) <- c("vargm", "list")
    return(vrgm)
    cat("\n")
}
