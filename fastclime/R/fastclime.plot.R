#-------------------------------------------------------------------------------#
# Package: fastclime                                                            #
# fastclime.generator(): graph visualization                                    #
# Authors: Haotian Pang, Han Liu and Robert Vanderbei                           #
# Emails: <hpang@princeton.edu>, <hanliu@princeton.edu> and <rvdb@princetonedu> #
# Date: April 22th 2016                                                           #
# Version: 1.4.1					                                            #
#-------------------------------------------------------------------------------#

fastclime.plot = function(G, epsflag = FALSE, graph.name = "default", cur.num = 1, location=NULL){
	gcinfo(FALSE)
	if(missing(location))	location = getwd()
	setwd(location)
        diag(G)=0
        Matrix(G,sparse=TRUE)
	g = graph.adjacency(as.matrix(G!=0), mode="undirected", diag=FALSE)
	layout.grid = layout.fruchterman.reingold(g)
	
   	if(epsflag == TRUE)	postscript(paste(paste(graph.name, cur.num, sep=""), "eps", sep="."), width = 8.0, height = 8.0)             
	par(mfrow = c(1,1))
	plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=2, vertex.label=NA)
	rm(g,location)	
   	gc()
   	if(epsflag == TRUE) dev.off()
}
