"pcol" <-
function(x, y = NULL, z = NULL, width = NULL, bins = 5, method = "pearson", permutations = 1000, alpha = 0.5, trace = FALSE, ...) {
	x <- as.matrix(x)
	N <- dim(x)[1]
	N <- N * (N-1) / 2 
    if(!is.null(z)){
    	z <- as.matrix(z)
		if(is.null(width)) {width <- ceiling(diff(range(z))/bins)}
		if(width == 1){
			z.nms <- unique(as.vector(z))
		}
		else{
			z <- floor(z/width)
			z.nms <- levels(as.factor(z*width))
		} 
    	if(is.null(y)){
    		strata <- unique(as.vector(z))
			stratan <- length(strata)
			tmp <- t(sapply(c(1:stratan), function(u) mantl(x, ifelse(z==strata[u], 1, 0), method=method, permutations=permutations, ...)[3:4]))
			nop <- summary(as.factor(z[row(z)>col(z)]))
			out <- data.frame(nop=nop, sapply(data.frame(tmp), unlist))
			rownames(out) <- z.nms
        	attr(out, "solo") <- TRUE
    	}
    	else {
			strata <- unique(as.vector(z))
			stratan <- length(strata)
			out <- data.frame(cbind(nop=1, statistic=1:stratan, signif=1))
			if (trace) {cat(stratan, "classes: ")}
			for (i in 1:stratan) {
				sub <- z==strata[i]
				tmp <- mantl(x, y, method=method, permutations=permutations, sub=sub, ...)
				out[i,3] <- as.numeric(tmp$statistic)
				out[i,2] <- as.numeric(tmp$signif)
				out[i,1] <- as.numeric(tmp$n)
				if (trace) {cat(paste(i,""))}
			}
    	}
    	alpha <- alpha/stratan
    	out$sig <- as.vector(symnum(out$signif, corr=FALSE, cutpoints = c(0,.001,.01,.05,.1,1), symbols = c("***","**","*",".","ns")))
    	rownames(out) <- strata
    	res <- list(call=match.call(), method=method, out=out, gesN=N, strata=stratan, permutations=permutations)
    	class(res) <- "pclist"
    	out <- res
    }
    else {
    	out <- mantl(x, y, method=method, permutations=permutations, ...)
    }
    out$call <- match.call()
    return(out)
}