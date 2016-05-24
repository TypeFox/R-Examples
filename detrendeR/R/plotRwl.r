plotRwl  =
function (Rwl, file.name="", biweight = TRUE, plot.mean=T, plot.spline=T, plot.y1=F, band.width=10, Spar=0.5, save.csv=F,...){
Year=as.numeric(rownames(Rwl))
TRW=Rwl
matplot(Year, TRW, ylab="Ring-width (mm)", type="l", lty=1, lwd=1, las=1, col="grey20",...)

if (plot.mean) {
 if (!biweight) 
        #    res = apply(x.ar, 1, mean, na.rm = TRUE)
             mean.rwl=apply(Rwl,1, mean, na.rm=T)
        else  mean.rwl=apply(Rwl,1, tbrm, C = 9)   #res = apply(x.ar, 1, tbrm, C = 9)

lines(Year, mean.rwl, col="red", lwd=3)}
if (plot.spline) {

SPLINE(as.vector(mean.rwl), bandwidth=band.width, p=Spar)->crn.spline
lines(Year, crn.spline, col="blue", lwd=2)}

if (plot.y1) abline(h=1)

std.dev<-apply(Rwl,1, sd, na.rm=T)
if(save.csv){
samp.depth = apply(TRW, 1, function(y) sum(!is.na(y)))
crn.raw<-data.frame(Year, mean.rwl, crn.spline, std.dev, samp.depth)
write.table(round(crn.raw,3), paste(file.name, "RAW.csv", sep=""), quote = FALSE, sep = ";", row.names =F)
}
}
