framework.stats <-
function (xx) 
{
    stats <- gx.stats(xx, display = FALSE, iftell = FALSE)
    table <- numeric(20)
    table[1] <- stats$stats[20]
    table[2] <- length(xx) - stats$stats[20]
    table[3:13] <- stats$stats[c(1, 3, 4, 5, 7, 10, 13, 15, 16, 
        17, 19)]
    table[14] <- stats$stats[27]
    table[15] <- stats$stats[28]
    table[16] <- stats$stats[21]
    table[17] <- stats$stats[22]
    table[18] <- stats$stats[23]
    table[19] <- stats$stats[25]
    table[14:19] <- signif(table[14:19], 4)
    table[20] <- stats$stats[26]
    invisible(table)
}
