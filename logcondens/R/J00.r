J00 <- function (x, y, v = 1){
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 0.005]
    z[II] <- z[II] * (exp(v * d[II]) - 1)/d[II]
    II <- (1:m)[abs(d) <= 0.005]
    z[II] <- z[II] *
           (v +
               d[II]*(v/2 +
                   d[II]*(v/6 +
                       d[II]*(v/24 + d[II]*v/120))))
    return(z)
}

# version prior to e-mails of LD of August 2010
# J00 <- function(x, y, v = 1){
#     m <- length(x)
#     z <- exp(x)
#     d <- y - x
#     II <- (1:m)[abs(d) > 10 ^ (-6)]
#     z[II] <- z[II] * (exp(v * d[II]) - 1) / d[II]
#     II <- (1:m)[abs(d) <= 10 ^ (-6)]
#     z[II] <- z[II] * (v + v ^ 2 * d[II]/2 + v ^ 3 * d[II] ^ 2 / 6)
#     return(z)
# }
