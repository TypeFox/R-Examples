"plot.scan.gwaa" <- 
function (x, y, ..., df=1, ystart=0, col=c("blue","green"), sort=TRUE, ylim, 
		delta = 1, main = getcall(x)) {
   df <- as.character(df)
   dfpars <- c("1","2","Pc1df","Pc2df")
   if (!(any(dfpars==df))) stop ("df parameter must be 1, 2, \"Pc1df\", or \"Pc2df\"")
   if (!is(x,"scan.gwaa")) stop("Plot must be done on an object of scan.gwaa-class")
   if (length(map(x)) != length(x[,"P1df"])) stop("length of map and scan points not equal!")
   #if (any(names(x) == "Pgw1df")) if (!is.null(x$Pgw1df)) {x$P1df <- x$Pgw1df; x$P2df <- x$Pgw2df}

   Pv <- x[,"P1df"]
   if (df=="2") {Pv <- x[,"P2df"];}
   else if (df=="Pc1df") {Pv <- x[,"Pc1df"];}
   else if (df=="Pc2df") {Pv <- x[,"Pc2df"];}

	if (sort) {
		newmap <- sortmap.internal(chromosome(x),map(x),delta=delta)
		mymap <- newmap$cummap
		mychromosome <- chromosome(x)[newmap$ix]
		Pv <- Pv[newmap$ix]
		newchnum <- newmap$chnum
		rm(newmap)
		gc()
	} else {
		mymap <- map(x)
		mychromosome <- chromosome(x)
	}

    Pv <- replace(Pv,(Pv<=0),1.e-16)
	maxy <- max(-log10(Pv),na.rm=TRUE)
	cargnams <- names(match.call())
	
	if (!missing(ylim)) {
		ylim <- ylim
	} else {
		ylim <- c(ystart,maxy)
	}

	if (dim(table(mychromosome))>1) {
		chind <- chrom.char2num(mychromosome) %% length(col)
		idxCH <- which(chind==0)
		plot(mymap[idxCH],-log10(Pv[idxCH]),main=main,xlab="Chromosome",
				ylab=expression(-log[10](P-value)), axes=FALSE, ylim=ylim, 
				xlim=c(min(mymap),max(mymap)), col=col[length(col)],...)
		for (colidx in c(1:(length(col)-1))) {
			idxCH <- which(chind==colidx)
			points(mymap[idxCH],-log10(Pv[idxCH]),col=col[colidx],...)
		}
		mxlog <- floor(max(-log10(Pv),na.rm=T))
		if (mxlog==0) mxlog <- max(-log10(Pv),na.rm=T)
		if (mxlog<1) 
			axis(2,at=ylim)
		else
			axis(2,at=c(ylim[1]:ylim[2]))
		chrs <- unique(mychromosome)
		chpos <- rep(NA,length(chrs))
		for (i in 1:length(chrs)) chpos[i] <- mean(mymap[mychromosome==chrs[i]])
		axis(1,at=chpos,labels=chrs)
	} else {
		plot(mymap,-log10(Pv),main=main,xlab="Map position",
				ylab=expression(-log[10](P-value)), 
				ylim=ylim, col=col[1], ...)
	}
}

