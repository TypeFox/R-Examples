
AppendCheckDigit.Verhoeff <- function(x) {

    x[is.na(x)] <- NA

    if (any(grepl('\\D', x))) {
        warning('Non-digit characters are disregarded in check digit calculation')
    }

    d <- matrix(
        c(
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            1, 2, 3, 4, 0, 6, 7, 8, 9, 5,
            2, 3, 4, 0, 1, 7, 8, 9, 5, 6,
            3, 4, 0, 1, 2, 8, 9, 5, 6, 7,
            4, 0, 1, 2, 3, 9, 5, 6, 7, 8,
            5, 9, 8, 7, 6, 0, 4, 3, 2, 1,
            6, 5, 9, 8, 7, 1, 0, 4, 3, 2,
            7, 6, 5, 9, 8, 2, 1, 0, 4, 3,
            8, 7, 6, 5, 9, 3, 2, 1, 0, 4,
            9, 8, 7, 6, 5, 4, 3, 2, 1, 0
        ),
        ncol=10,
        byrow=TRUE
    )

    p <- matrix(
        c(
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            1, 5, 7, 6, 2, 8, 3, 0, 9, 4,
            5, 8, 0, 3, 7, 9, 6, 1, 4, 2,
            8, 9, 1, 6, 0, 4, 3, 5, 2, 7,
            9, 4, 5, 3, 1, 2, 6, 8, 7, 0,
            4, 2, 8, 6, 5, 7, 3, 9, 0, 1,
            2, 7, 9, 3, 8, 0, 6, 4, 1, 5,
            7, 0, 4, 6, 9, 1, 3, 2, 5, 8
        ),
        ncol=10,
        byrow=TRUE
    )

    inv <- matrix(
        c(0, 4, 3, 2, 1, 5, 6, 7, 8, 9),
        ncol=10,
        byrow=TRUE
    )

    y <- sapply(x, function(x) {
        y <- gsub('\\D', '', x)
        if (y %in% c('', NA, NaN)) {
            return(x)
        } else {
            c <- 0
            n <- as.integer(rev(unlist(strsplit(y, ''))))
            for (i in seq(n)) {
                c <- d[c+1, p[i %% 8 + 1, n[i]+1]+1]
            }
            return(paste(x, inv[c+1], sep=''))
        }
    })

    return(y)
}

VerifyCheckDigit.Verhoeff <- function(x) {

    x[is.na(x)] <- NA

    if (any(grepl('\\D', x))) {
        warning('Non-digit characters are disregarded in check digit calculation')
    }

    d <- matrix(
        c(
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            1, 2, 3, 4, 0, 6, 7, 8, 9, 5,
            2, 3, 4, 0, 1, 7, 8, 9, 5, 6,
            3, 4, 0, 1, 2, 8, 9, 5, 6, 7,
            4, 0, 1, 2, 3, 9, 5, 6, 7, 8,
            5, 9, 8, 7, 6, 0, 4, 3, 2, 1,
            6, 5, 9, 8, 7, 1, 0, 4, 3, 2,
            7, 6, 5, 9, 8, 2, 1, 0, 4, 3,
            8, 7, 6, 5, 9, 3, 2, 1, 0, 4,
            9, 8, 7, 6, 5, 4, 3, 2, 1, 0
        ),
        ncol=10,
        byrow=TRUE
    )

    p <- matrix(
        c(
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            1, 5, 7, 6, 2, 8, 3, 0, 9, 4,
            5, 8, 0, 3, 7, 9, 6, 1, 4, 2,
            8, 9, 1, 6, 0, 4, 3, 5, 2, 7,
            9, 4, 5, 3, 1, 2, 6, 8, 7, 0,
            4, 2, 8, 6, 5, 7, 3, 9, 0, 1,
            2, 7, 9, 3, 8, 0, 6, 4, 1, 5,
            7, 0, 4, 6, 9, 1, 3, 2, 5, 8
        ),
        ncol=10,
        byrow=TRUE
    )

    y <- sapply(x, function(x) {
        y <- gsub('\\D', '', x)
        if (y %in% c('', NA, NaN)) {
            return(FALSE)
        } else {
            c <- 0
            n <- as.integer(rev(unlist(strsplit(y, ''))))
            for (i in seq(n)) {
                c <- d[c+1, p[(i-1) %% 8 + 1, n[i]+1]+1]
            }
            return(c == 0)
        }
    })

    return(y)
}
