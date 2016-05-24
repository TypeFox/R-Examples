J20 <- function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 0.02]
    z[II] <- 2 * z[II] *
        (exp(d[II]) - 1 - d[II] - d[II]^2/2)/(d[II]^3)
    II <- (1:m)[abs(d) <= 0.02]
    z[II] <- z[II] *
        (1/3 + 
            d[II]*(1/12 +
                d[II]*(1/60 +
                    d[II]*(1/360 + d[II]/2520))))
    return(z)
}

# version prior to e-mails of LD of August 2010
# J20 <- function(x, y){
#     m <- length(x)
#     z <- exp(x)
#     d <- y - x
#     II <- (1:m)[abs(d) > 10^(-6)]
#     z[II] <- 2 * z[II] * (exp(d[II]) - 1 - d[II] - d[II] ^ 2 / 2)/(d[II] ^ 3)
#     II <- (1:m)[abs(d) <= 10 ^ (-6)]
#     z[II] <- z[II] * (1 / 3 + d[II] / 12 + d[II] ^ 2 / 60 + d[II] ^ 3 / 360)
#     return(z)
# }
