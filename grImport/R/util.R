
# This is very single-byte-encoding-specific, but that
# should be ok because that is all I am expecting to
# get out of tracing a PostScript file
tidyString <- function(string) {
    # Convert any non-printing char to a full stop
    rawString <- charToRaw(string)
    nonPrinting <-  rawString < 32
    if (any(nonPrinting)) {
        rawString[nonPrinting] <- as.raw(46)
        rawToChar(rawString)
    } else {
        string
    }
}
# Test:
# tidyString(rawToChar(as.raw(1:40)))

readLTY <- function(ltyString) {
    tc <- textConnection(ltyString)
    lty <- scan(tc, quiet=TRUE)
    close(tc)
    lty
}

# Go from numeric vector for lty to string
fixLTY <- function(lty, lwd) {
    if (length(lty)) {
        # Minimum allowed is 1
        paste(as.hexmode(pmax(1, round(lty/lwd))),
              collapse="")
    } else {
        "solid"
    }
}
