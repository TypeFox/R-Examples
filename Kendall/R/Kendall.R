"Kendall" <-
function(x, y)
{
    #calculate Kendall rank correlation
    if (length(x) != length(y)) stop("length of inputs muct be equal!")
    if (length(x) < 3) stop("length(x)<3")
    tau <- 0.0
    ptau <- 0.0
    sltau <- 0.0
    score <- 0.0
    varscore <- 0.0
    denom <- 0.0
    i <- !(is.na(x)|is.na(y))
    x <- x[i]
    y <- y[i]
    iws <- numeric(length(x))
    ifault <- 0.0
    outF <- .Fortran("tauk2",
        as.single(x),
        as.single(y),
        as.integer(length(x)),
        as.single(tau),
        as.single(ptau),
        as.single(sltau),
        as.single(score),
        as.single(varscore),
        as.single(denom),
        as.integer(iws),
        as.integer(ifault),
        PACKAGE="Kendall")
    tau <- outF[[4]]
    sl <- outF[[6]]
    sc <- outF[[7]]
    var.sc <- outF[[8]]
    denom <- outF[[9]]
    ier <- outF[[11]]
    if(ier != 0.) {
        cat(paste("WARNING: Error exit, tauk2. IFAULT = ", ier), fill = TRUE)
    }
    ans <- list(tau = tau, sl = sl, S = sc, D = denom, varS = var.sc)
    oldClass(ans) <- "Kendall"
    ans
}

