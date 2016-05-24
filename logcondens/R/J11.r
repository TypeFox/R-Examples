J11 <- function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 0.02]
    z[II] <- z[II] *
        (d[II]*(exp(d[II]) + 1) - 2*(exp(d[II]) - 1))/(d[II]^3)
    II <- (1:m)[abs(d) <= 0.02]
    z[II] <- z[II] *
        (1/6 + 
            d[II]*(1/12 +
                d[II]*(1/40 +
                    d[II]*(1/180 + d[II]/1008))))
    return(z)
}


# version prior to e-mails of LD of August 2010
# J11 <- function(x, y){
#     m <- length(x)
#     z <- exp(x)
#     d <- y - x
#     II <- (1:m)[abs(d) > 10^(-6)]
#     z[II] <- z[II] * (exp(d[II]) * (d[II] - 2) + 2 + d[II])/(d[II] ^ 3)
#     II <- (1:m)[abs(d) <= 10^(-6)]
#     z[II] <- z[II] * (1 / 6 + d[II] / 12 + d[II] ^ 2 / 60 + d[II] ^ 3 / 180)
#     return(z)
# }
