gx.summary <-
function (xx, log = log, iftell = iftell) 
{
    nxx <- length(xx)
    table <- numeric(15)
    stats <- gx.stats(xx, display = FALSE, iftell = iftell)
    table[1] <- stats$stats[20]
    table[2] <- nxx - stats$stats[20]
    table[3:7] <- stats$stats[c(1, 7, 10, 13, 19)]
    table[8] <- signif(stats$stats[21], 4)
    table[9] <- signif(stats$stats[22], 4)
    if (log) {
        table[10] <- stats$stats[29]
        table[11] <- stats$stats[31]
        table[12] <- stats$stats[32]
        table[13] <- stats$stats[31]/sqrt(stats$stats[20])
        ci <- qt(0.975, stats$stats[20] - 1) * table[13]
        table[14] <- 10^(table[10] - ci)
        table[15] <- 10^(table[10] + ci)
        table[10] <- 10^table[10]
    }
    else {
        table[10] <- stats$stats[23]
        table[11] <- stats$stats[25]
        table[12] <- stats$stats[26]
        table[13] <- stats$stats[25]/sqrt(stats$stats[20])
        ci <- qt(0.975, stats$stats[20] - 1) * table[13]
        table[14] <- table[10] - ci
        table[15] <- table[10] + ci
    }
    table[3:15] <- signif(table[3:15], 4)
    invisible(table)
}
