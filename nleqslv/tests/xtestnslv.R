
library(nleqslv)

# function to replace small number with OK or if not with NZ
# this is to avoid differences in the Fnorm column between machines/cpu/os/compilers

# the test is for checking that testnslv (still) works as expected

fixsmall <- function(x) {
    z <- ifelse(x < .Machine$double.eps^(2/3), "OK","NZ")
    z <- ifelse(is.na(z), "NA", z)
    z
}

dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}
xstart <- c(0.5,0.5)
fstart <- dslnex(xstart)
z <- testnslv(xstart,dslnex)
zfn <- z$out[,"Fnorm"]
z$out[,"Fnorm"] <- fixsmall(zfn)
z

# this will encounter an error
xstart <- c(2.0,0.5)
fstart <- dslnex(xstart)
z <- testnslv(xstart,dslnex)
zfn <- z$out[,"Fnorm"]
z$out[,"Fnorm"] <- fixsmall(zfn)
z
