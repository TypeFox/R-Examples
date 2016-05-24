
pamr.plotfdr <- function(fdrfit,  call.win.metafile=FALSE){

##if(call.win.metafile){win.metafile()}

om=fdrfit$results[,"Number of significant genes"]>0

na.min=function(x){min(x[!is.na(x)])}
na.max=function(x){max(x[!is.na(x)])}
plot(fdrfit$results[om,"Number of significant genes"],fdrfit$results[om,"Median FDR"],log="x",
xlab="Number of genes called significant",
ylab="False discovery rate (median and 90th percentile)",type="b",  
ylim=c(na.min(fdrfit$results[om,"Median FDR"]), na.max(fdrfit$results[om,"90th percentile of FDR"])))
x=fdrfit$results[om,"Number of significant genes"]
xlim <- range(x)
barw <- abs((log(x)))*1.2
upper=fdrfit$results[om,"90th percentile of FDR"]
lower=fdrfit$results[om,"Median FDR"]
segments(x, upper, x, lower, lty=2)
segments(x - barw, upper, x + barw, upper, lty=2)

axis(3,at=fdrfit$results[om,"Number of significant genes"], labels=round(fdrfit$results[om,"Threshold"],2))

mtext("Threshold", 3, 2, cex = 1.0)

if (call.win.metafile) {
    savePlot("", type="wmf")
}
dev.off()

return()
}
