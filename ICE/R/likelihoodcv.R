"likelihoodcv" <-
function (p, data, m=101) 
{
    inner <- inmost(data)
    cvsamples <- function(data, innermost) {
        n.inner <- length(innermost[, 1])
        datacv <- vector("list", length = n.inner)
        for (i in 1:n.inner) {
            datacv[[i]] <- data[(data[, 1] != innermost[i, 1]) & 
                (data[, 2] != innermost[i, 2]), ]
        }
        list(datacv = datacv, n.inner=n.inner)
    }
    iccdf <- function(xvalues, data, h, m) {
        yvalues <- NULL
        f.ic <- ickde(data, h, m = m)
        for (x in xvalues) {
            if (x < f.ic$x[1]) {
                yvalues <- c(yvalues, 0)
            }
            else {
                if (x > f.ic$x[m]) {
                  yvalues <- c(yvalues, 1)
                }
                else {
                  delta <- f.ic$x[2] - f.ic$x[1]
                  yvalues <- c(yvalues, (sum(f.ic$y[f.ic$x <= 
                    x]) * delta))
                }
            }
        }
        yvalues
    }
    samples <- cvsamples(data, inner)
    n.inner <- samples$n.inner
    samples <- samples$datacv
    ecdf <- numeric(n.inner)
    for (j in 1:n.inner) {
        ecdf[j] <- diff(iccdf(inner[j, 1:2], samples[[j]], h = p, m=m))
    }
    -sum(log(ecdf[abs(ecdf)>.00000001]))
}

