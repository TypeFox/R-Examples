setGeneric("hrBootstrap", function(x, rep=100, plot=TRUE, level=95, levelMax=100, unin='km', unout='m2', ...){standardGeneric("hrBootstrap")})
setMethod("hrBootstrap", 
	  signature=c(x="SpatialPoints"),
	  definition=function(x, ...){
		  if (class(x)=="SpatialPoints") nlocs <- length(x) else nlocs <- n.locs(x)
		  j <- round(exp(log(10)*seq(log10(5),log10(nlocs),by=0.1)))
		  message(paste(j, collapse=" "))
		  crds<-coordinates(x)
		  dimnames(crds)<-NULL
		  locs <- lapply(j, function(js,coords) {replicate(n=rep,SpatialPoints(coords[sample(1:nlocs,size=js,replace=T),]) )}, coords=crds)
		  if (requireNamespace("adehabitatHR", quietly = TRUE)) {
			  mcps <- lapply(locs, lapply, adehabitatHR::mcp.area, percent=level, plotit=F,unin=unin,unout=unout)
			  quant <- lapply(lapply(mcps, unlist), quantile, probs=seq(0,1,.25))
			  quantLines <- as.data.frame(do.call("rbind", quant))
			  rownames(quantLines) <- j
			  qll <- do.call("list", quantLines)
			  a <- list(2,3,1,3,2)
			  col <- list(5,2,1,2,5)
			  b <- list(qll,a,col)

			  if (class(x)=="SpatialPoints"){
				  hline <- as.numeric(adehabitatHR::mcp.area(x, percent=levelMax, plotit=F, unin=unin, unout=unout),lty=7)
			  } else {
				  hline <- as.numeric(adehabitatHR::mcp.area(move2ade(x), percent=levelMax, plotit=F, unin=unin, unout=unout),lty=7)
			  }         
		  } else {
			  stop("The package adehabitatHR can't be loaded, probably it needs to be installed")
		  }
		  if (plot){
			  plot(x=j, quantLines[,ncol(quantLines)], xlab="number of coordinates for mcp", ylab=paste("habitat area (",unout,")"), type="n", ylim=range(c(quantLines,hline)), ...)
			  abline(h=hline)
			  lapply(1:ncol(quantLines), function(i, j, b) lines(x=j, y=b[[1]][[i]],lty=b[[2]][[i]],col=b[[3]][[i]],lwd=2), j=j, b=b)
		  }
		  return(quantLines)
	  })


setMethod("hrBootstrap", 
	  signature=c(x=".MoveTrackStack"),
	  definition=function(x, ...){
		  return(lapply(split(x), hrBootstrap, rep=rep, plot=plot, level=level, levelMax=levelMax, unin=unin, unout=unout, ...))
	  })

