treetext<-function (hc,pval, col = c(2, 3), print.num = TRUE, float = 0.01, 
    cex = NULL, font = NULL) 
{
    axes <- hc2axes(hc)
    usr <- par()$usr
    wid <- usr[4] - usr[3]
    bp <- as.character(round(pval * 100))
    rn <- as.character(1:length(pval))
    bp[length(bp)] <- "pvalue*100"
    rn[length(rn)] <- "edge #"
    a <- text(x = axes[, 1], y = axes[, 2] + float * wid, bp, 
        col = col[1], pos = 2, offset = 0.3, cex = cex, font = font)
    if(print.num){
        a <- text(x = axes[, 1], y = axes[, 2], rn, col = col[2], 
            pos = 4, offset = 0.3, cex = cex, font = font)
	}
}
