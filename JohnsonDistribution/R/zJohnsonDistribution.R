"zJohnsonDistribution" <-
function(s, ITYPE, GAMMA, DELTA, XLAM, XI)
{
#Johnson-to-normal transformation
   stopifnot(length(ITYPE)==1, length(GAMMA)==1, length(DELTA)==1, length(XLAM)==1, length(XI)==1)
   ifault <- 0.0
   SNV <- 0.0
   n <- length(s)
   z <- numeric(n)
   for (i in 1:n) {
    AJV = s[i]
    outF <- .Fortran("SUBSNV",
        as.single(SNV),
        as.single(AJV),
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
    z[i] <- outF[[1]]
    }
    z
}
