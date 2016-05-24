wtd.boxplot.stats <-
function(x, weights=NULL, coef = 1.5, do.conf=TRUE, do.out=TRUE)
{

    nna <- !is.na(x)
    n <- sum(nna)                       # including +/- Inf
#   stats <- stats::fivenum(x, weights=weights, na.rm = TRUE) # is the new call
# the previous lines needs to be uncommented fot inclusion in the R distribution
# and the next line has to be deleted
    stats <- wtd.fivenum(x, weights=weights, na.rm = TRUE) # is the call for the test version
    iqr <- diff(stats[c(2, 4)])
    if(coef < 0) stop("'coef' must not be negative")
    if(coef == 0)
do.out <- FALSE
    else {                              # coef > 0
out <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
if(any(out[nna])) stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
    conf <- NULL
    if(do.conf && is.null(weights)) conf <- stats[3] + c(-1.58, 1.58) * iqr / sqrt(n)
    if(do.conf&& !is.null(weights))
    conf <- stats[3] + c(-1.58, 1.58) * iqr * sqrt(sum((weights/sum(weights))^2))
    nn<-ifelse(is.null(weights),n,sum(weights))
    list(stats = stats, n = nn, conf = conf,
 out = if(do.out) x[out & nna] else numeric(0))
}

