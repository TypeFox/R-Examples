summary.criterionRkh <-
function(object,...){
x <- object
if (!inherits(x, "criterionRkh")) stop("use only with \"criterionRkh\" objects")

cat("Trace criterion estimated by boostrap", "\n")
cat(" \n")

cat(paste("Reduction method performed:", x$method),"\n")
cat(paste("Number of observations:", x$n),"\n")

cat(" \n")
cat("Value of K tried:")
cat(" \n")
cat(x$K)
cat(" \n")
cat("Value of H tried:")
cat(" \n")
cat(x$H)
cat(" \n")
cat(" \n")

cat("Result of the bootsrap estimate of the trace criterion R(k,h) \n" )

tmp <- x$Rkhbootmean
row.names(tmp) <- x$H
colnames(tmp) <- x$K
cat("\n")
prmatrix(signif(tmp,3))
cat("\n")
}

