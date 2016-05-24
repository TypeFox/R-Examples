pitHist <-
function( fit, ensembleData, dates=NULL)
{

 if (class(fit) == "ensembleBMAgamma0" || class(fit) == "fitBMAgamma0") {
   randomizeATzero <- TRUE
 }
 else {
   randomizeATzero <- FALSE
 }

 if (!is.null(dates)) {
   d <- match(dates, ensembleValidDates(ensembleData))
   PIT <- pit( fit[dates], ensembleData[d,], randomizeATzero = randomizeATzero)
 }
 else {
   PIT <- pit( fit, ensembleData, randomizeATzero = randomizeATzero)
 }

 k <- ensembleSize(ensembleData)

 hist(PIT, breaks = (0:(k+1))/(k+1),  prob=TRUE, xlab="", xaxt="n",
     ylab = "", main = "Probability Integral Transform")
 axis(1, at = seq(0, to = 1, length=11), labels = (0:10)/10)
 abline(h = 1, lty = 2)
 invisible(PIT)
}

