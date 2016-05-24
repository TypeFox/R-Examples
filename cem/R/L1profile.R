# ref.cut : number of cutpoints of other variables when cutting on one variable
# max.cut : maximal number of cut points
# min.cut : minimal number of cut points

L1.profile <- function(group, data, drop = NULL, min.cut = 2, max.cut = 12, 
weights, plot = TRUE, add = FALSE, col = "red", 
lty = 1, M=100, useCP=NULL, grouping=NULL, progress=TRUE) 
{
    if (min.cut < 2) min.cut <- 2
    min.cut <- as.integer(min.cut)
    max.cut <- as.integer(max.cut)
	if(max.cut<min.cut) max.cut <- min.cut+1
    if (!is.null(drop)) {
        drop <- unique(drop)
        dropped <- match(drop, colnames(data))
        dropped <- dropped[!is.na(dropped)]
        if (length(dropped) > 0) 
		data <- data[-dropped]
    }
    vnames <- colnames(data)
    nv <- length(vnames)
    idx <- sapply(1:nv, function(x) is.numeric(data[[x]]) | is.integer(data[[x]]))
    numvar <- vnames[which(idx)]
    catvar <- vnames[which(!idx)]
    R <- sapply(numvar, function(x) range(data[x], na.rm = TRUE), 
				simplify = FALSE)
    br <- vector(length(numvar), mode = "list")
    names(br) <- numvar
    out <- NULL
	if(is.null(grouping))
	 grouping <- "NONE"
	outCP <- vector(0, mode="list")
	outGR <- vector(0, mode="list")
	for (i in catvar) br[i] <- NULL
	ns <- 0
	cat("\n")
	if(!is.null(useCP)){
		nCP <- length(useCP)
		if(progress){pb <- txtProgressBar(min = 1, max = nCP, initial = 1, style = 3)}
#		s <- round(nCP*.1)
		for(k in 1:nCP){
			if(progress){setTxtProgressBar(pb, k)}
			cp <- useCP[[k]]
			ncp <- length(cp)
			for(j in 1:ncp)
				br[[ names(cp)[j] ]] <-  cp[[j]]
			tmp <- L1.meas(group = group, data = data, breaks = br, 
						   weights = weights, grouping=grouping)$L1
			names(tmp) <- names(useCP)[k]
			out <- c(out, tmp)
			outCP[[length(outCP)+1]] <-  br
			outGR[[length(outGR)+1]] <-  grouping
#			if(k %% s == 0){
#				ns <- ns+10
#				cat(sprintf("[%2d%%]",ns))
#			} else {
#				cat(".")
#			}			
		}
		if(progress){close(pb)}	
	} else {
#		s <- round(M*.1)
		if(progress){pb <- txtProgressBar(min = 1, max = M, initial = 1, style = 3)}
 	 	for(m in 1:M){
			if(progress){setTxtProgressBar(pb, m)}
			theta <- sample(min.cut:max.cut, length(numvar),replace=TRUE)
			names(theta) <- numvar
			for (i in numvar) 
				br[[i]] <- unique(seq(R[[i]][1], R[[i]][2], length = theta[i]))
			tmp <- L1.meas(group = group, data = data, breaks = br, 
						   weights = weights, grouping=grouping)$L1
			names(tmp) <- sprintf("rnd(%d)",m)
			out <- c(out, tmp)
			outCP[[length(outCP)+1]] <-  br
			outGR[[length(outGR)+1]] <-  grouping
#			if(m %% s == 0){
#				ns <- ns+10
#				cat(sprintf("[%2d%%]",ns))
#			} else {
#				cat(".")
#			}
			
		}
		if(progress){close(pb)}
	}
	
	names(outCP) <- names(out)
	names(outGR) <- names(out)
    idx <- order(out)
	if(is.null(useCP)){
	 out <- out[idx]
	 outCP <- outCP[idx]
	 outGR <- outGR[idx]
	}
	medianL1 <- median(out)
	medianCP <- outCP[[ which(sort(out) >= medianL1)[1] ]]
	medianGR <- outGR[[ which(sort(out) >= medianL1)[1] ]]
	
    val <- list(L1 =out, CP = outCP, GR=outGR, medianL1=medianL1, medianCP=medianCP, medianGR=medianGR)
    class(val) <- "L1profile"
	cat("\n")
    if (plot) 
	plot(val, lty = lty, col = col, add = add)
    return(invisible(val))
}


plot.L1profile <- function (x, add = FALSE, cex.axis = 0.3, ...) 
{
    if (!add) {
        plot(unclass(x$L1), type = "n", axes = F, xlab = "cut points", 
			 ylab = expression(L[1]), ylim = c(0, 1), ...)
        axis(1, 1:length(x$L1), names(x$L1), padj = 0.5, las = 3, xpd = TRUE, 
			 cex.axis = cex.axis)
        axis(2)
    }
    lines(x$L1, ...)
}
