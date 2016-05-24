gx.stats <-
function (xx, xlab = deparse(substitute(xx)), display = TRUE, iftell = TRUE) 
{
    rgv <- numeric(13)
    probs <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 
        0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.98, 0.99, 1)
    qtls <- numeric(19)
    stats <- numeric(32)
    temp.x <- remove.na(xx, iftell = iftell)
    if (temp.x$n < 2) {
        rgv[1] <- temp.x$n
        rgv[2:13] <- qtls[1:19] <- NA  
        if (temp.x$n == 1) {
            rgv[9] <- rgv[8] <- rgv[4] <- qtls[c(1,10,19)] <- temp.x$x[1]
            if(temp.x$x[1] > 0) rgv[10] <- log10(temp.x$x[1])
        }
    }
    else {    
        x <- sort(temp.x$x[1:temp.x$n])
        qtls <- quantile(x, probs = probs)
        rgv[1] <- temp.x$n
        rgv[2] <- mad(x)
        rgv[3] <- (qtls[13] - qtls[7]) * 0.7413
        rgv[4] <- mean(x)
        rgv[5] <- var(x)
        rgv[6] <- sqrt(rgv[5])
        rgv[7] <- round((100 * rgv[6])/rgv[4], 2)
        j <- qbinom(0.025, temp.x$n, 0.5)
        if (j < 1) j <- 1
        rgv[8] <- x[j]
        rgv[9] <- x[temp.x$n - j + 1]
        if (qtls[1] >= 0) {
            logx <- log10(x)
            rgv[10] <- mean(logx)
            rgv[11] <- var(logx)
            rgv[12] <- sqrt(rgv[11])
            rgv[13] <- round((100 * rgv[12])/rgv[10], 2)
        }
        else rgv[10:13] <- NA
    }
    if (display) {
        cat("\n Summary Statistics Display for:", xlab, "\n\n Data Set N =  ", 
            rgv[1], "\n Minimum =     ", qtls[1], "\t\tMaximum =", 
            qtls[19], "\n Median =      ", signif(qtls[10], 4), "\t\tMAD Est =", 
            signif(rgv[2], 4), "\n\t\t\t\tIQR Est =", signif(rgv[3], 
                4), "\n 95% CI for the Median =\t  ", rgv[8], 
            "to", rgv[9], "\n\n Mean =        ", signif(rgv[4], 
                4), "\t\tS.D. =   ", signif(rgv[6], 4), "\n Variance =    ", 
            signif(rgv[5], 4), "\t\tC.V. % = ", rgv[7], "\n")
        cat("\n Maximum Value         ", qtls[19], "\n", "99th Percentile       ", 
            qtls[18], "\n", "98th Percentile       ", qtls[17], 
            "\n", "95th Percentile       ", qtls[16], "\n", "90th Percentile       ", 
            qtls[15], "\n", "80th Percentile       ", qtls[14], 
            "\n", "3rd Quartile (75th)   ", qtls[13], "\n", "70th Percentile       ", 
            qtls[12], "\n", "60th Percentile       ", qtls[11], 
            "\n", "Median (50th)         ", qtls[10], "\n", "40th Percentile       ", 
            qtls[9], "\n", "30th Percentile       ", qtls[8], 
            "\n", "1st Quartile (25th)   ", qtls[7], "\n", "20th Percentile       ", 
            qtls[6], "\n", "10th Percentile       ", qtls[5], 
            "\n", " 5th Percentile       ", qtls[4], "\n", " 2nd Percentile       ", 
            qtls[3], "\n", " 1st Percentile       ", qtls[2], 
            "\n", "Minimum Value         ", qtls[1], "\n\n")
    }
    stats[1:19] <- signif(qtls[1:19], 4)
    stats[20:32] <- rgv[1:13]
    invisible(list(stats = stats))
}
