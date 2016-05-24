"linlogtics" <- function (x, mu, alpha) 
{
    maxorigx <- max(x)
    minorigx <- min(x)
    ticsl <- c(-alpha)
    tics <- c(-alpha)
    cntmin <- -alpha
    
    while (cntmin > minorigx) {
        cntmin <- cntmin * 10
        ticsl <- append(ticsl, cntmin)
        tics <- append(tics, -alpha - (alpha * log10(-cntmin/alpha)))
    }
    ticsl <- append(sort(ticsl), c(0, alpha))
    tics <- append(sort(tics), c(0, alpha))
    cntmax <- alpha
    while (cntmax < maxorigx) {
        cntmax <- cntmax * 10
        ticsl <- append(ticsl, cntmax)
        tics <- append(tics, alpha + (alpha * log10(cntmax/alpha)))
    }
    ## new x values as column 1
    ## new x labels as colum 2
    ret<-cbind(tics, ticsl)

    ret
}
