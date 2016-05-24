"FitJohnsonDistribution" <-
function(XBAR, SD, RB1, BB2)
{
    #Fit Johnson Distribution
    stopifnot(length(XBAR)==1, length(SD)==1, length(RB1)==1, length(BB2)==1)
    ITYPE <- 0
    GAMMA <- 0.0
    DELTA <- 0.0
    XLAM <- 0.0
    XI <- 0.0
   ifault <- 0.0
    outF <- .Fortran("JNSN",
        as.single(XBAR),
        as.single(SD),
        as.single(RB1),
        as.single(BB2),
        as.integer(ITYPE),
        as.single(GAMMA),
        as.single(DELTA),
        as.single(XLAM),
        as.single(XI),
        as.integer(ifault),
        PACKAGE="JohnsonDistribution")
    ITYPE <- outF[[5]]
    GAMMA <- outF[[6]]
    DELTA <- outF[[7]]
    XLAM <- outF[[8]]
    XI <- outF[[9]]
    ier <- outF[[10]]
    if(ier != 0.)
        cat(paste("WARNING: Error exit, JNSN. IFAULT = ", ier), fill = TRUE)
    ans <- c(ITYPE, GAMMA, DELTA, XLAM, XI, ier)
    names(ans)<-c("ITYPE", "GAMMA", "DELTA", "XLAM", "XI", "IFAULT")
    ans
}

