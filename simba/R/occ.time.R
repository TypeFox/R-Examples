"occ.time" <- 
function(x, y, times=NULL, adjust=TRUE, gen.occ=FALSE, perc=TRUE, nc.acc=FALSE, ...) {	
	if(is.null(times)){
	   bac <- occ.tmp(x, y, adjust = adjust, gen.occ=gen.occ, perc=perc, nc.acc=nc.acc)
	}
	else{
	   dat <- x
	   strata <- unique(times)
	   stratan <- length(strata)
	   id.nms <- outer(strata, strata, paste, sep=".")
	   id.nms <- id.nms[row(id.nms) < col(id.nms)]
	   out <- lapply(c(1:stratan), function(x) dat[(times==strata[x]),])
	   out <- lapply(out, mama)
	   id <- outer(c(1:stratan), c(1:stratan), paste)
       id <- strsplit(id[row(id) < col(id)], " ")
       id.l <- c(1:length(id))
       id <- cbind(as.numeric(sapply(id, function(x) x[1])), as.numeric(sapply(id, function(x) x[2])))
       outk <- lapply(id.l, function(x) occ.tmp(out[[(id[x,1])]], out[[(id[x,2])]], adjust=adjust, gen.occ=gen.occ, perc=perc, nc.acc=nc.acc))
       out <- lapply(outk, function(x) x[[1]])
       n.plots <- sapply(outk, function(x) x[[2]])
       n.spec <- sapply(outk, function(x) x[[3]])
       n.spec1 <- sapply(outk, function(x) x[[4]])
       n.spec2 <- sapply(outk, function(x) x[[5]])
       n.spec1o <- sapply(outk, function(x) x[[6]])
       n.spec2o <- sapply(outk, function(x) x[[7]])
       spec.nms1o <- lapply(outk, function(x) x[[8]])
       spec.nms2o <- lapply(outk, function(x) x[[9]])
       stats <- cbind(n.plots, n.spec, n.spec1, n.spec2, n.spec1o, n.spec2o)
       out.kat <- lapply(out, attr, "names")
       max.t <- max(sapply(out.kat, function(x) max(as.numeric(x))))
       min.t <- min(sapply(out.kat, function(x) min(as.numeric(x))))
       bac <- matrix(c(min.t:max.t), (max.t-min.t)+1, 1)
       bac <- data.frame(bac)
       out <- lapply(out, data.frame)
       out <- lapply(out, function(x) merge(data.frame(bac), data.frame(x), by.x="bac", by.y=0, all=TRUE))
       bac <- data.frame(bac, sapply(out, function(x) x[ ,2]))
       rownames(bac) <- bac[,1]
       bac <- bac[,-1]
       names(bac) <- id.nms
       bac <- t(bac)
       stats <- data.frame(stats)
       rownames(stats) <- rownames(bac)
	    }
	   res <- list(bac=bac, stats=stats, spec.nms1o=spec.nms1o, spec.nms2o=spec.nms2o)
	   class(res) <- "occtmp"
	return(res)
}