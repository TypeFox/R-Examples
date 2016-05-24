rowVars <-
function(x, ...) {
 n <- NCOL(x)
 newn <- n-rowSums(is.na(x))
 val <-  (rowSums(x^2, ...) - newn * rowMeans(x, ...)^2)/(newn - 1)
 val <- replace(val,which((newn)<2),NA)
 val
}

