signpred <- function(matbin,pred.lablength=max(sapply(rownames(matbin),nchar)),labsize=1,plotsize = 12){
if(is.null(rownames(matbin))){rownames(matbin) <- paste("x",1:nrow(matbin),sep="")}
lll=max(sapply(rownames(matbin),nchar))
text = "no"
ncol <- ncol(matbin)
nrow <- nrow(matbin)
plotsize = plotsize/2.54
mcol = max(matbin)
if (ncol > nrow) {
        wx <- plotsize
        wy <- (plotsize)/ncol * nrow
}
else {
        wy <- plotsize
        wx <- (plotsize)/nrow * ncol
}
m.colsize = max(strwidth(colnames(matbin), units = "inches"))
m.rowsize = max(strwidth(rownames(matbin), units = "inches"))
cellsize = wx/ncol
if (substr(text, 1, 1) == "i") 
        s <- as.character(max(matbin))
else s = "A"
lettersize = strwidth(s, units = "inches")
clratio = cellsize/lettersize
mm <- max(m.colsize, m.rowsize)
op <- par(las=3,mar = c(2, 2, 1, 1) + 0.1, mgp = c(2, 1, 0))
on.exit(par(op))
visweb(t(matbin),type="None",labsize=labsize,square="b",box.col="grey25",prednames=FALSE,clear=FALSE)
text(.5+0:(length(rownames(matbin))-1),-lll+1,rownames(matbin),cex=0.4 * labsize *clratio,srt=90)
}
