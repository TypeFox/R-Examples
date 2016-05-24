'plot.angledistlist' <- function(x, makepdf=FALSE,...){

		if(makepdf)pdf("All distributions.pdf",onefile=TRUE)
		dists <- names(x$allfits)
		for(d in dists){
				plot.angledist(x$allfits[[d]],...)
		}
		if(makepdf)dev.off()
}