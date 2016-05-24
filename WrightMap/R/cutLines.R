cutLines <- function(cutpoints = NULL,cut.left = 0, cut.right = 1, cut.lab.text = NULL, cut.lab.adj = c(0,1),...) {
	
	if(!is.null(cutpoints)) {
		
		width <- par("usr")[2]
	
	segments(cut.left * width,cutpoints,cut.right * width, col = 'grey70', lwd = 1.5, lend = 2)
		
		text(width *(.01  + cut.lab.adj[1] * .98),cutpoints + (.5 - cut.lab.adj[2])*.1,labels = cut.lab.text,cex= .6, adj = cut.lab.adj)
		
		}
}