timeSeq.genesetscreeplot = function(timeSeq.geneset){
  graphics.off()
  
  geneset = timeSeq.geneset
  n = length(geneset$NPDE.ratio)  
  NPDE.order = order(-geneset$NPDE.ratio)  
  NPDE.list = data.frame(geneset.name = geneset$geneset.names[NPDE.order], 
	  				     ratio = geneset$NPDE.ratio[NPDE.order])  
  plot(x = 1 : n, y = NPDE.list$ratio, type = "b", pch = 21, col = "red", xaxt = "n", 
	   lty = 2, main = "Ratios of Genesets", xlab = "Geneset", ylab = "Ratio")    
  axis(1, 1 : n, NPDE.list$geneset.name[1 : n], col.axis = "blue")
  out = NPDE.list
  out
}
