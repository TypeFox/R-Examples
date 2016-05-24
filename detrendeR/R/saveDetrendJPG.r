saveDetrendJPG = function (rwl, detrend, folderName = "Detrend", work.dir=NULL, detrend.method="", select.series =1:(ncol(rwl)) ){

if (is.null(work.dir)) work.dir = getwd()
dir.create(folderName, showWarnings = FALSE)
setwd(folderName)

for (i in select.series){
seriesnames=colnames(rwl)
yr.vec = as.numeric(rownames(rwl))
jpeg(paste(seriesnames[i], ".jpg", sep=""),width = 1200, height = 600, quality=100)
plot(yr.vec, rwl[,i], type="l", xlab = "Years",ylab = "", main = seriesnames[i], las = 1, col="blue")
mtext(paste(detrend.method, sep=""), line=0.5, side =3, adj =1, cex=0.90, col="blue",font=1)
mtext("Detrender", line=0.5, side =3, adj =0, cex=0.90, col="blue",font=1)
lines(yr.vec, detrend[,i], col=2)
dev.off()
}
setwd(work.dir)
}