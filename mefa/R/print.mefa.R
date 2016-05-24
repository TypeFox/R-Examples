`print.mefa` <-
function(x, nlist=10, ...)
{
#cat("\nCall: ", deparse(x$call), "\n", sep="")
cat("\nAn object of class 'mefa' containing\n\n", sep="")
indiv <- if (sum(x$xtab) == 1) " individual" else " individuals"
cat(" $ xtab: ", sum(x$xtab), indiv, " of ", ncol(x$xtab), " taxa in ", nrow(x$xtab), " samples,\n", sep="")
cat(" $ segm: ", sep="")

if (!is.null(x$segm)) {
    cat(length(x$segm), sep="")
    if (attr(x, "nested")) cat(" nested segments:\n", sep="") else cat(" (non-nested) segments:\n", sep="")
    dimn <- if (dim(x)[3] > nlist) {
        c(paste(dimnames(x)[[3]][1:nlist], rep(", ", nlist), sep=""), "...")
        } else {
        paste(dimnames(x)[[3]], c(rep(", ", dim(x)[3]-1), ""), sep="")}
    cat("         ", dimn, ",\n", sep="")
    } else cat("1 (all inclusive) segment,\n", sep="")

prov <- if (is.null(x$samp)) "not provided,\n" else {
    paste("provided (", NCOL(x$samp), " variables),\n", sep="")}
cat(" $ samp: table for samples ", prov, sep="")
prov <- if (is.null(x$taxa)) "not provided." else {
    paste("provided (", NCOL(x$taxa), " variables).", sep="")}
cat(" $ taxa: table for taxa ", prov, sep="")
cat("\n\n")
invisible(x)
}

