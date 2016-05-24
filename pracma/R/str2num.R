##
##  s t r 2 n u m . R
##


str2num <- function(S) {
    s1 <- strTrim(S)
    ls <- nchar(s1)
    if (substr(s1, ls, ls) == ';') {
        s1 <- sub(';$', '', s1)
        prit <- FALSE
    } else {
        prit <- TRUE
    }

    s1 <- sub('^\\[', '', s1)
    s1 <- sub('\\]$', '', s1)
    s1 <- gsub(',', ' ', s1)

    s2 <- strsplit(s1, ';')[[1]]

    m  <- length(s2)
    L1 <- scan(text=s2[1], quiet = TRUE)
    n  <- length(L1)

    if (m > 1) {
        for (i in 2:m) {
            Li <- scan(text=s2[i], quiet = TRUE)
            if (n != length(Li))
                stop("All rows in Argument 's' must have the same length.")
            L1 <- rbind(L1, Li)
        }
    }
    L2 <- unname(L1)
    if (any(is.na(L2)) || isempty(L2)) L2 <- c()

    if (prit) print(L2)
    invisible(L2)
}


num2str <- function(A, fmt = 3) {
    stopifnot(is.numeric(A), length(fmt) == 1)
    if (is.numeric(fmt))
        fmt = paste("%.", round(fmt), "f", sep = '')

    dm <- dim(A)
    a1 <- sprintf(fmt, A)       # a2 <- as.numeric(a1)
    if (!is.null(dm)) {
        dim(a1) <- dm           # dim(a2) <- dm
    }
    return(a1)
}
