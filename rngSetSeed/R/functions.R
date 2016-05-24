setVectorSeed <- function(vseed)
{
    RNGkind("Mersenne-Twister")
    stopifnot(length(.Random.seed) == 626)
    stopifnot(.Random.seed[2] == 624)
    .Random.seed[3:626] <<- generateInitialization(vseed, 624)
    invisible(NULL)
}

generateInitialization <- function(vseed, m)
{
    if (any(vseed != floor(vseed))) stop("Vector seed should have integer components")
    if (any(vseed < 0 | vseed >= 2^32)) stop("Vector seed should have components in [0, 2^32-1]")
	vseed <- c(vseed, length(vseed))
    s <- numeric(8*ceiling(length(vseed)/8))
    s[seq.int(along.with=vseed)] <- vseed
    m4 <- 4*ceiling(m/4)
    .C("getVectorSeed",
        length(s),
        as.double(s),
        as.integer(m4),
        out=integer(m4),
        PACKAGE="rngSetSeed")$out[seq.int(length.out=m)]
}

