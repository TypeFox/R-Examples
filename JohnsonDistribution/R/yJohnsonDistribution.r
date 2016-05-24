"yJohnsonDistribution" <-
function(z, ITYPE, GAMMA, DELTA, XLAM, XI)
{
#Normal to Johnson transformation
   stopifnot(length(ITYPE)==1, length(GAMMA)==1, length(DELTA)==1,
    length(XLAM)==1, length(XI)==1)
   ifault <- 0.0
   AJV <- 0.0
   n <- length(z)
   s <- numeric(n)
   for (i in 1:n) {
    zi <- z[i]
    outF <- .Fortran("SUBAJV",
        as.single(AJV),
        as.single(zi),
        as.integer(ITYPE),
        as.single(GAMMA),
        as.single(DELTA),
        as.single(XLAM),
        as.single(XI),
        as.integer(ifault),
        PACKAGE="JohnsonDistribution")
    ier <- outF[[8]]
    if(ier != 0.)
        cat(paste("WARNING: Error exit, SUBAJV. IFAULT = ", ier), fill = TRUE)
    s[i] <- outF[[1]]
    }
    s
}
