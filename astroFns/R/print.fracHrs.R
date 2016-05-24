# Print method for fracHrs class
# AH 2012.05.20-2012.06.15

print.fracHrs <- function(x, digits = NULL, quote = TRUE, na.print = NULL,
                          print.gap = NULL, right = FALSE, max = NULL,
                          useSource = TRUE, ...) {

    if (!inherits(x, 'fracHrs')) stop('*** Wrong class for print.fracHrs ***')

    # Set up format with number of decimal places in time
    # Use option('digit.secs') for this to match e.g. POSIXct()
    pl <- getOption('digits.secs')
    if (is.null(pl)) pl <- 0
    pl <- max(pl, 0)
    pl <- min(pl, 6)
    # Output format
    if (pl == 0) {
        fmt <- '%02d:%02d:%02d'
    } else {
        pl <- 3+pl + (pl)/10
        fmt <- paste('%02d:%02d:%0#', pl, 'f', sep='')
    }

    # Output: convert to string with hours, minutes, seconds
    txt <- character(length(x))
    for (i in 1:length(x)) {
        h <- unclass(x[i])
        m <- (abs(h)%%1) * 60
        s <- (m%%1) * 60

        h <- trunc(h)
        m <- trunc(m)
        s <- round(s, pl)
        # Convert 60 secs to 0 secs with incremented minutes, etc.
        if(s >= 60) {
            s <- s-60
            m <- m+1
        }
        if(m >= 60) {
            m <- m-60
            h <- h+1
        }
        if(h >= 24) {
            h <- h-24
            x <- x+1
        }
        txt[i] <- sprintf(fmt, h, m, s)
    }

    cat(txt, '\n', sep='  ', fill=TRUE)

}

