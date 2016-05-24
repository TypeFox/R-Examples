plot.vargm <- function (x, smooth = FALSE, bdw = NULL, 
                        follow.time = NULL, points = TRUE, ...) 
{
    vargm <- x
    svar <- vargm$svar
    sigma2 <- vargm$sigma2
    if (is.null(follow.time)) {
        left <- min(svar[, 1])
        right <- max(svar[, 1])
    }
    else {
        left <- follow.time[1]
        right <- follow.time[2]
    }
    if (is.null(bdw)) {
        bdw <- (right - left)/10
    }
    nbdw <- (right - left)/bdw
    Mid <- c()
    Mean <- c()
    Count <- c()
    for (it in 1:(nbdw)) {
        lt <- left + (it - 1) * bdw
        rt <- left + it * bdw
        mid <- lt + (rt - lt)/2
        if (length(svar[svar[, 1] < rt & svar[, 1] >= lt, 2]) > 
            0) {
            mean <- mean(svar[svar[, 1] < rt & svar[, 1] >= lt, 
                2])
            Count <- c(Count, length(svar[svar[, 1] < rt & svar[, 
                1] >= lt, 2]))
        }
        else {
            mean <- NA
            Count <- c(Count, 0)
        }
        Mid <- c(Mid, mid)
        Mean <- c(Mean, mean)
    }
    mean.line <- cbind(Mid[!is.na(Mean)], Mean[!is.na(Mean)])
    Count <- Count[!is.na(Mean)]
    if (points == T) {
        plot(svar[, 1], svar[, 2], xlab = "u", ylab = "Variogram", 
            pch = ".", ...)
        if (smooth == T) {
            lines(smooth.spline(svar[, 1], svar[, 2]), lty = 1, 
                lwd = 1.3)
        }
        else {
            lines(mean.line, lwd = 1.3, lty = 1)
        }
        abline(h = sigma2, lwd = 1.3, lty = 1)
    }
    else {
        if (smooth == T) {
            plot(smooth.spline(svar[, 1], svar[, 2]), , xlab = "u", 
                ylab = "Variogram", ...)
        }
        else {
            plot(mean.line, type = "l", xlab = "u", ylab = "Variogram", 
                ...)
        }
        abline(h = sigma2, lwd = 2, lty = 1)
        abline(h = min(mean.line[mean.line[, 1] <= right & 
            mean.line[, 1] >= left, 2]), lty = 2)
        abline(h = max(mean.line[mean.line[, 1] <= right & 
            mean.line[, 1] >= left, 2]), lty = 2)
    }
    cat("\n")
}
