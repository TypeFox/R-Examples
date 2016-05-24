print.BDtest <-
function(x, ...)
{

cat("Confidence intervals for binary diagnostic tests.\n")

cat(attr(x$INDAT, which="caption"), "\n")
print.data.frame(x$INDAT,...)

cat(attr(x$SESPDAT, which="caption"), "\n")
print.data.frame(x$SESPDAT,...)

cat(attr(x$PPVNPVDAT, which="caption"), "\n")
print.data.frame(x$PPVNPVDAT,...)

}

