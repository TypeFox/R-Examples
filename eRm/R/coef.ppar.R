`coef.ppar` <-
function(object, extrapolated = TRUE, ...) {
   x <- object$theta.table[,1]
   if(!extrapolated) x[object$theta.table[,3]] <- NA
   names(x) <- rownames(object$theta.table)
   x
}
