##
##   a g m . R  Arithmetic-geometric Mean
##


agm <- function(a, b) {
    stopifnot(is.numeric(a) || is.complex(a),
              is.numeric(b) || is.complex(b))
    if (length(a) > 1 || length(b) > 1)
        stop("Arguments 'a', 'b' must be real or complex scalars.")
    if (is.numeric(a) && a < 0 || is.numeric(b) && b < 0) {
        a <- as.complex(a)
        b <- as.complex(b)
    }

    while (a != b) {
        a1 <- (a + b) / 2
        b1 <- sqrt(a * b)
        a <- a1
        b <- b1
    }
    a
}

