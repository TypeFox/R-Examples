plot.CCM <-
function(x, y, index, no.plot = FALSE, ...) {
   KK = split(x[index,], y)
   if (no.plot) return(KK)
   boxplot(KK, ...)
}

