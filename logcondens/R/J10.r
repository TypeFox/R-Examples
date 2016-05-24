J10 <- function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 0.01]
    z[II] <- z[II] * (exp(d[II]) - 1 - d[II])/(d[II]^2)
    II <- (1:m)[abs(d) <= 0.01]
    z[II] <- z[II] *
        (1/2 + 
            d[II]*(1/6 +
                d[II]*(1/24 +
                    d[II]*(1/120 + d[II]/720))))
    return(z)
}

# version prior to e-mails of LD of August 2010
# J10 <- function(x, y){
#     m <- length(x)
#     z <- exp(x)
#     d <- y - x
#     II <- (1:m)[abs(d) > 10 ^ (-6)]
#     z[II] <- z[II] * (exp(d[II]) - 1 - d[II]) / (d[II] ^ 2)
#     II <- (1:m)[abs(d) <= 10 ^ (-6)]
#     z[II] <- z[II] * (1 / 2 + d[II]/6 + d[II] ^ 2 / 24)
#     return(z)
# }
