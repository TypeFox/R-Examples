timeSeq.screeplot = function(timeSeq.obj, type = c("barplot", "lines")) {
	#n = min(n, length(obj$NPDE.ratio))
	graphics.off()
    obj_sorted = timeSeq.obj$sorted
    obj_sorted = obj_sorted$NPDE_list
    gene.length = sum(!is.na(obj_sorted[,2]))
    gene.index = !(is.na(obj_sorted[,2]))
    if (type == "barplot") {
        obj_sorted$genenames = factor(obj_sorted$genenames, 
                                      levels = unique( as.character(obj_sorted$genenames)))    
        barchart(ratio ~ genenames, data = obj_sorted, 
                 main = "Ratios of Genes", xlab = "Gene", ylab = "Ratio",
                 col = "aliceblue")    
    } else if (type == "lines") {
        plot(x = 1 : gene.length, y = obj_sorted$ratio[gene.index],
             type = "b", pch = 21, col = "red", xaxt = "n", lty = 2, 
             main = "Ratios of Genes",  xlab = "Gene", ylab = "Ratio")    
        axis(1, 1 : gene.length, obj_sorted$genenames[gene.index], col.axis = "blue")
	} else {
		cat("The type of plot must be bar plot or line graph.\n")
		return(cat("ERROR!"))
	} 
}
