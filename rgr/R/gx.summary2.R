gx.summary2 <-
function (xx, xname = deparse(substitute(xx))) 
{
    nxx <- length(xx)
    table <- numeric(28)
    stats <- gx.stats(xx, display = FALSE, iftell = FALSE)
    table[1] <- stats$stats[20]
    table[2] <- nxx - stats$stats[20]
    table[3] <- stats$stats[23]
    temp <- qt(0.975, stats$stats[20] - 1)/sqrt(stats$stats[20])
    ci <- temp * stats$stats[25]
    table[4] <- table[3] - ci
    table[5] <- table[3] + ci
    table[6] <- stats$stats[25]
    table[7] <- stats$stats[26]
    table[8] <- stats$stats[29]
    table[9] <- stats$stats[31]
    ci <- temp * stats$stats[31]
    table[10] <- 10^stats$stats[29]
    table[11] <- 10^(stats$stats[29] - ci)
    table[12] <- 10^(stats$stats[29] + ci)
    table[13] <- stats$stats[21]
    table[14] <- stats$stats[22]
    table[15] <- stats$stats[10]
    table[16] <- stats$stats[27]
    table[17] <- stats$stats[28]
    table[18:28] <- stats$stats[c(1, 3:5, 7, 10, 13, 15:17, 19)]
    table[3:6] <- signif(table[3:6], 4)
    table[8:28] <- signif(table[8:28], 4)
    cat("\n  Extended Summary Stats for: ", xname, "\n  N and number of NAs:\t\t", 
        table[1], " & ", table[2], "\n  Arithmetic Mean and 95% CLs:\t", 
        table[3], " & ", table[4], " <-> ", table[5], "\n  SD and CV%:\t\t\t", 
        table[6], " & ", table[7], "%\n  Geometric Mean and 95% CLs:\t", 
        table[10], " & ", table[11], " <-> ", table[12], "\n  Log10 Mean and SD:\t\t", 
        table[8], " & ", table[9], "\n  Median and 95% CLs:\t\t", 
        table[15], " & ", table[16], " <-> ", table[17], "\n  MAD and IQR estimates of SD:\t", 
        table[13], " & ", table[14], "\n  Percentiles: Min  2  5  10  25 - 50 - 75  90  95  98  Max\n\t", 
        table[18], " ", table[19], " ", table[20], " ", table[21], 
        " ", table[22], " - ", table[23], " - ", table[24], " ", 
        table[25], " ", table[26], " ", table[27], " ", table[28], 
        "\n\n", sep = "")
    invisible(table)
}
