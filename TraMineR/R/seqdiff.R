###########################
## Locate the difference in sequence between groups
###########################
seqdiff <- function(seqdata, group, cmprange=c(0, 1),
					seqdist_arg=list(method="LCS", norm=TRUE),
					with.missing=FALSE, weighted=TRUE, squared=FALSE) {

	if (!inherits(seqdata, "stslist")) {
			stop("seqdata should be a stslist, see seqdef")
	}
	seqdist_arg$with.missing <- with.missing
	slenE <- ncol(seqdata)
	ret <- NULL
	startAt <- 1
	#Range where we compare
	totrange <- max(startAt, 1-cmprange[1]):min(slenE, slenE-cmprange[2])
	ret <- list()
	name_column <- ""
	if (inherits(group, "stslist")) {
		name_column <- alphabet(group)
	} else {
		group=factor(group)
		name_column <- levels(group)
	}
	num_column <- length(name_column)
	ret$stat <- matrix(NA, nrow=length(totrange), ncol=5)
	rownames(ret$stat) <- colnames(seqdata)[totrange]
	colnames(ret$stat) <- c("Pseudo F", "Pseudo Fbf", "Pseudo R2", "Bartlett", "Levene")
	ret$discrepancy <- matrix(0, nrow=length(totrange), ncol=(num_column+1))
	rownames(ret$discrepancy) <- colnames(seqdata)[totrange]
	colnames(ret$discrepancy) <- c(name_column, "Total")
	
	weights <- attr(seqdata, "weights")
	if (!weighted) {
		weights <- NULL
	}

	for (i in totrange) {
		gc()
		srange=c((i+cmprange[1]):(i+cmprange[2]))
		if (inherits(group, "stslist")) {
			cmpbase <- group[, i]
			cmpbase[cmpbase==attr(group, "void") |cmpbase==attr(group, "nr")] <- NA
		} else {
			cmpbase <- group
		}
		#Getting complete case on range and group var
		subseq <- seqdata[, srange]
		if (!with.missing) {
			subseq2 <- subseq
			subseq2[subseq==attr(seqdata, "void") |subseq==attr(seqdata, "nr")] <- NA
			seqok <- complete.cases(cmpbase, subseq2)
		}
		else {
			seqok <- complete.cases(cmpbase)
		}
		seqdist_arg$seqdata <- subseq[seqok, ]
		## Computing distance on range
		sdist <- suppressMessages(do.call(seqdist, args=seqdist_arg))
		tmp <- dissassoc(sdist, cmpbase[seqok], R=0, weights=weights[seqok], weight.permutation="diss", squared=squared)
		ret$stat[i,rownames(tmp$stat)] <- tmp$stat[,1]
		ret$discrepancy[i,rownames(tmp$groups)] <- tmp$groups$discrepancy
	}
	attr(ret, "xtstep") <- attr(seqdata, "xtstep")
	class(ret) <- "seqdiff"
	return(ret)
}

###########################
## Print method for seqdiff
###########################
print.seqdiff <- function(x, ...) {
	message("\nStatistics:")
	print(x$stat, ...)
	message("\nDiscrepancies:")
	print(x$discrepancy, ...)
}

###########################
## Plot method for seqdiff
###########################
plot.seqdiff <- function(x, stat="Pseudo R2", type="l", ylab=stat, xlab="",
							legendposition="top", ylim=NULL, xaxt=TRUE, col=NULL, xtstep=NULL, ...) {
	if(is.null(xtstep)){
		xtstep <- ifelse(!is.null(attr(x, "xtstep")), attr(x, "xtstep"), 1)
	}
	if(length(stat)==1){
    	if (stat %in% c("Variance", "discrepancy", "Residuals","residuals")) {
    		nbstates=ncol(x$discrepancy)
    		if(is.null(col)) {
    			if (nbstates <= 8) cpal <- brewer.pal(nbstates, "Accent")
    			else if (nbstates > 8 & nbstates <= 12) cpal <- brewer.pal(nbstates, "Set3")
    		} else {
    			cpal <- col
    		}
    		if (stat %in% c("Residuals", "residuals")) {
    			toplot <- x$discrepancy*(1-x$stat[, "Pseudo R2"])
    		} else {
    			toplot <- x$discrepancy
    		}
    		if (is.null(ylim)) {
    			ylim=c(min(toplot), max(toplot))
    		}
    		plot(1:nrow(x$discrepancy), x$discrepancy[, ncol(x$discrepancy)], type=type, ylab=ylab, xlab=xlab, xaxt="n", col=cpal[nbstates], ylim=ylim, ...)

    		for (i in 1:(ncol(x$discrepancy)-1)) {
    			lines(toplot[, i], type=type, col=cpal[i], ...)
    		}
    		legend(legendposition, fill = cpal, legend = colnames(x$discrepancy))
			if(xaxt){
				tpos <- seq(1,nrow(x$discrepancy), xtstep)
				axis(1, at=tpos-0.5, labels=rownames(x$stat)[tpos])
			}
    		return(invisible())
    	}
	
	    ##if(length(stat)==1){
		if(is.null(col)){
			col <- c("black")
		}
		if (stat %in% colnames(x$stat)) {
			plot(x$stat[, stat], type=type, ylab=ylab, xlab=xlab, xaxt="n",col=col, ...)
		}
		else {
			stop(" [!] 'stat' argument should be one of ", paste(c("discrepancy", colnames(x$stat)), sep=", "), ".")
		}
		
	}else if (length(stat)==2) {
		if(sum(stat %in% colnames(x$stat)) < 2 ){
			stop(" [!] The two values of the 'stat' argument should be one of ", paste(colnames(x$stat), sep=", "), ".")
		}
        for (i in 1:2){
            if (!(stat[i] %in% colnames(x$stat))){
          			
            }
        }
		if(is.null(col)){
			col <- c("red", "blue")
		}
		plot(x$stat[, stat[1]], type=type, ylab=ylab, xlab=xlab, xaxt="n",col=col[1],axes=FALSE, ...)
		axis(2,col=col[1],col.axis=col[1])
		oldpar <- par(new=TRUE)
		on.exit(par(oldpar))
		plot(x$stat[, stat[2]],type=type,xlab="",ylab="", xaxt="n",col=col[2],axes=FALSE, ...)
		axis(4,col=col[2],col.axis=col[2])
		legend(legendposition, fill = col, legend =stat)
	}
	else{
		stop(" [!] Too many values for the 'stat' argument (max 2)")
	}
	if(xaxt){
		tpos <- seq(1,nrow(x$discrepancy), xtstep)
		axis(1, at=tpos-0.5, labels=rownames(x$stat)[tpos])
	}
}
