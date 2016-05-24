library(lattice)



panel.meansdplot <- function(x,y,mean.pch=15,mean.cex=1,mean.col="blue",...){
	meansdplot.stats <- function(x,...){
	    nna <- !is.na(x)
	    n <- sum(nna)
	    s <- sd(x, na.rm=TRUE)
	    m <- mean(x, na.rm=TRUE)
	    stats <- c(NA,m-s,m,m+s,NA)
	    list(stats=stats, n=n, NULL, c())
	}

	panel.bwplot(x=x,y=y,pch=mean.pch,cex=mean.cex,col=mean.col,box.width=0,stats=meansdplot.stats,do.out=FALSE,...);
}

panel.violinm <- function(x,y,horizontal=TRUE, grid=TRUE, mean.col="blue",violin.col="transparent", ...){
	if( grid ){
		if( horizontal ){
			panel.grid(h=0, v=-1)
		}else{
			panel.grid(h=-1, v=0)
		}
	}
	panel.violin(x=x,y=y,col=violin.col,horizontal=horizontal,...)
	panel.meansdplot(x=x,y=y,mean.col=mean.col,horizontal=horizontal,...)
}

violinmplot <- function(x, data, ...){
	bwplot( x=x, data=data, panel=panel.violinm, ...)
}
