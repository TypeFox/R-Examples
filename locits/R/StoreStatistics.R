StoreStatistics <-
function () 
{
ans <- .C("StoreStatistics", lfound = as.double(0), lstored=as.double(0), loutside=as.double(0), PACKAGE="locits")

loutside <- ans$loutside
lstored <- ans$lstored
lfound <- ans$lfound

cat("Number calculated outside framework: ", loutside, "\n")
cat("Number calculated then stored: ", lstored, "\n")
cat("Number found in store: ", lfound, "\n")

cat("\n")

cat("Overall % calculated: ", 100*(loutside+lstored)/(lfound+loutside+lstored), "\n")
cat("% outside framework: ", 100*loutside/(loutside+lfound+lstored), "\n")
}
