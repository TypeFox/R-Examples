
seqclustname <- function(seqdata, group, diss, weighted=TRUE, perc=FALSE){
	
	if(weighted){
		ww <- attr(seqdata, "weights")
	} else{
		ww <- NULL
	}
	xx <- disscenter(diss, group=group, medoids.index="first", allcenter = FALSE, weights=ww)
	
	suppressMessages(ll <- seqformat(seqdata[xx, ], from="STS", to="SPS", compressed=TRUE, SPS.out=list(xfix="", sdsep="/")))
	if(perc){
		if(is.null(ww)){
			tt <- xtabs(~group)
		}
		else{
			tt <- xtabs(ww~group)
		}
		tt <- tt[names(xx)] 
		tt <- round(prop.table(tt)*100, 1)
		ll <- paste(ll, " (",tt, "%)", sep="")
	}
	return(factor(group, levels=names(xx), labels=ll))
}
