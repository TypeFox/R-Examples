##
##  g a m m a z . R
##


gammaz <- function(z) {
    if (!is.numeric(z) && !is.complex(z))
        stop("Argument 'z' must be a numeric or complex vector.")

    z  <- c(z); zz <- z
    f  <- complex(length(z))

    p  <- which(Re(z) < 0)
    if (length(p) > 0) z[p] <- -z[p]

    # Lanczos approximation
    g  <- 607/128

    cc <- c(  0.99999999999999709182,     57.156235665862923517,
            -59.597960355475491248,       14.136097974741747174,
             -0.49191381609762019978,      0.33994649984811888699e-4,
              0.46523628927048575665e-4,  -0.98374475304879564677e-4,
              0.15808870322491248884e-3,  -0.21026444172410488319e-3,
              0.21743961811521264320e-3,  -0.16431810653676389022e-3,
              0.84418223983852743293e-4,  -0.26190838401581408670e-4,
              0.36899182659531622704e-5)

    z <- z - 1
    zh <- z + 0.5
    zgh <- zh + g
    # trick for avoiding FP overflow above z=141
    zp <- zgh^(zh*0.5)

    ss <- 0.0
    for (k in (length(cc)-1):1)
        ss <- ss + cc[k+1] / (z + k)

    # sqrt(2Pi)
    sq2pi <- 2.5066282746310005024157652848110;
    f <- (sq2pi * (cc[1]+ss)) * ((zp * exp(-zgh)) * zp);

    f[z == 0 | z == 1] <- 1.0

    # adjust for negative real parts
    if (length(p) > 0)
        f[p] <- -pi/(zz[p] * f[p] * sin(pi*zz[p]))

    # adjust for negative poles
    p <- which(round(zz) == zz & Im(zz) == 0 & Re(zz) <= 0)
    if (length(p) > 0)
        f[p] <- Inf

    return(f)
}


# lgammaz <- function(z) {
#     if (!is.numeric(z) && !is.complex(z))
#         stop("Argument 'z' must be a numeric or complex vector.")
# 
#     dimz <- dim(z)
#     z  <- c(z) + 0i
#     zz <- z
#     f  <- complex(length(z))
# 
#     p  <- which(Re(z) < 0)
#     if (length(p) > 0) z[p] <- -z[p]
# 
#     # Lanczos approximation
#     g  <- 607/128
# 
#     cc <- c(  0.99999999999999709182,     57.156235665862923517,
#             -59.597960355475491248,       14.136097974741747174,
#              -0.49191381609762019978,      0.33994649984811888699e-4,
#               0.46523628927048575665e-4,  -0.98374475304879564677e-4,
#               0.15808870322491248884e-3,  -0.21026444172410488319e-3,
#               0.21743961811521264320e-3,  -0.16431810653676389022e-3,
#               0.84418223983852743293e-4,  -0.26190838401581408670e-4,
#               0.36899182659531622704e-5)
# 
#     s <- 0
#     for (k in length(cc):2)
#         s <- s + cc[k]/(z + (k-2))
# 
#     zg <- z + g - 0.5
#     s2pi <- 0.9189385332046727417803297
#     f <- (s2pi + log(cc[1]+s)) - zg + (z-0.5) * log(zg)
# 
#     f[z == 1 | z == 2] <- 0.0
# 
#     if (length(p) > 0) {
#         lpi <- 1.14472988584940017414342735 + pi*1i
#         f[p] <- lpi - log(zz[p]) - f[p] - log(sin(pi*zz[p]))
#     }
# 
#     p <- which(round(zz) == zz & Im(zz) == 0 & Re(zz) <= 0)
#     if (length(p) > 0)
#         f[p] <- Inf
# 
#     dim(f) <- dimz
#     return(f)
# }

