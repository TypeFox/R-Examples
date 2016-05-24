e.scal<-function(x, k = 1, tc = NULL)
{
    scal <- function(k, x)
    {
        xx <- exp(k * x)
        return(xx / rowSums(xx))
    }
    scal2 <- function(k, x, tc)
    {
        xx <- scal(k, x)
        ec <- factor(max.col(xx), levels = seq(along = colnames(xx)), 
            labels = colnames(xx))
        return(mean(ec != tc))
    }
    if(!is.null(tc)) {
        k <- optimize(scal2, c(0, 1000), x = x, tc = tc)$minimum
    }
    x <- scal(k, x)
    return(list(sv = x, k = k))
}
