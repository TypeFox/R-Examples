spcdisc <- function(x,sort=FALSE)
{
   shannon <- function(y) {
       y <- as.numeric(y)
       frac <- y/sum(y)
       comp <- sum((-1 * frac * log(frac))[frac>0])
       comp <- 1 - (comp/log(length(y)))
       comp
   }
   tmp <- apply(x,1,shannon)
   if (sort) tmp <- tmp[rev(order(tmp))]
   tmp
}
