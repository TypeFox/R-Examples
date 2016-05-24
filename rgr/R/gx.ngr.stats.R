gx.ngr.stats <-
function (xx) 
{
    table <- numeric(29)
    nxx <- length(xx)
    stats <- gx.stats(xx, display = FALSE, iftell = FALSE)
    table[1] <- stats$stats[20]                    # N
    table[2] <- nxx - stats$stats[20]              # NAs
    table[3] <- stats$stats[23]                    # Mean
    table[4] <- stats$stats[25]                    # SD
    table[5] <- gx.ngr.skew(xx)                    # Skew
    table[6] <- stats$stats[26]                    # CV%
    table[7] <- 10^stats$stats[29]                 # Geom Mean
    table[8] <- stats$stats[10]                    # Median
    table[9] <- stats$stats[21]                    # MAD
    table[10] <- round(100 * table[9]/table[8], 2) # Robust CV%
    table[11:29] <- stats$stats[1:19]              # Table <- stats
    table[3:29] <- signif(table[3:29], 4)          # Table <- %iles
    invisible(table)
}
