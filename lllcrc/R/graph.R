#' Stack local capture pattern frequencies for plotting
#' 
#' Used in the plot functions
#' 
#' @param ydens A matrix or data.frame with as many columns as there are
#' capture patterns.  Each row is the relative frequency (possibly non-integer)
#' of each pattern.
#' @author Zach Kurtz
#' @export stackydens
stackydens = function(ydens){ 
  # location = code/graphics.R
  ystack = ydens
  for(k in 2:ncol(ydens)){ystack[, k] = ystack[, k-1] + ystack[, k]}
  return(ystack)}


#' Plot LLLMs
#' 
#' Generates stacked relative frequency diagrams for local log-linear or VGAM
#' capture recapture models
#' 
#' The capture pattern relative frequencies are plotted in a stacked form,
#' summing to 1 over each vertical cross-section.  On top, the rate of
#' missingness is plotted for the all-zero capture pattern, along with 95
#' percent confidence curves if a \code{boots} entry is included in \code{x},
#' the first argument
#' 
#' @aliases plot.lllcrc
#' @param x Output of the function for LLLMs, \code{lllcrc}, or the CRC wrapper
#' function for VGAM, \code{vgam.crc}
#' @param cont.var The name (in the form \code{x.dis ...} -- see
#' \code{format.data}) of a continuous variable that will be used for the
#' x-axis of the plot.  By default, it is \code{x.con.1}, the first continuous
#' variable
#' @param selection See the same argument in \code{extract.CI}
#' @param main Plot title
#' @param xlab x-axis label
#' @param label.offset Controls the distance between the capture-pattern labels
#' and corresponding curves
#' @param text.size Controls text size of capture pattern labels
#' @param x.range The range of the primary continuous variable that is included
#' in the plot
#' @param subtitle An optional subtitle
#' @param padj.adj Controls the vertical positioning of the subtitle if it is
#' defined
#' @param lr Has value 1 or -1; controls the left-right (or right-left)
#' alternating pattern of the capture pattern labels
#' @param lr.global NULL by default; if specified, it must be -1 or 1, causing
#' all capture pattern labels to appear on the left or the right, without
#' alternating
#' @param ylim Optional argument of the form c(a,b), where a and b are numbers.
#' @param ...  Additional parameters to be passed into \code{plot}
#' @return Returns nothing, but makes a plot if you're paying attention
#' @author Zach Kurtz
#' @rdname plot
#' @method plot lllcrc
#' @export
plot.lllcrc = function(x, cont.var = "x.con.1", 
		selection = NULL,
                main = "Set the main argument for the title",
                xlab = "x",
                label.offset = 0.01,
                text.size = 0.7,
                x.range = NULL, 
		subtitle = NULL,
                padj.adj = -0.3,
		lr = 1,
		lr.global = NULL,
		ylim = NULL, ...
	)
{
	xax = x$dat[,cont.var]
 	y = x$hpi
	# order y by average capture pattern frequency
	ynam = colnames(y); y = data.frame(y); names(y) = ynam
	y = y[, rev(order(colSums(y)))]
	Y = names(y)
 	ny = length(Y)
	stk = stackydens(y)
	stk$x = xax
	stk$pi0 = x$dat$pi0
	top = 1+max(stk$pi0)
	mct = x$dat$mct
	if(!is.null(x$boots)){
		stk[, c("l", "u")] = t(apply(x$boots$loc.est, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))/mct
		top = min(1+max(stk$u), 1+2*max(stk$pi0))
	}
	if(!is.null(selection)) {
		if(is.numeric(selection)) {
			stk = stk[selection,]
		}else{
			keepers = which(rowSums(x$dat[,selection, drop = FALSE]) == length(selection))
			stk = stk[keepers,]
		}
	}
	stk = stk[order(stk$x),]
	nx = nrow(stk)
	if(is.null(ylim)) ylim = c(0, top + 0.1)
	if(is.null(x.range)) x.range = diff(range(stk$x))
	# initialize plot
	plot(stk$x, stk[,1], , bty = "n", ylim = ylim,
		xlim = c(min(stk$x)-0.1*x.range, max(stk$x)+0.1*x.range),
		type = "n", col = 1, ylab = "Stacked relative frequency", 
		xlab = xlab, main = main, ...)
	# put in the curves
	for(k in 1:ny) lines(stk$x, stk[,k])
	# put labels on each curve
	for(k in ny:1){
		lr = -1*lr; if(!is.null(lr.global)) lr = lr.global
		if(lr > 0){
		this.x = max(stk$x) + label.offset*(x.range)
		text(Y[k], x = this.x, y = stk[nx,k], cex = text.size, adj = 0)
		}else{
		this.x= min(stk$x) - label.offset*(x.range)
		text(Y[k], x = this.x, y = stk[1,k], cex = text.size, adj = 1)
		}
	}#end for
	# add imputation
	lines(stk$x, stk$pi0+1, lwd = 2) # put in the final conditional curve
	k = nchar(names(stk)[1])
	text(paste(rep("0",k), collapse = ""), x = min(stk$x) - label.offset*(x.range), 
		y = 1 + stk$pi0[1], cex = text.size, adj = 1)
	if(!is.null(subtitle)) mtext(subtitle, padj=padj.adj, cex=text.size*1.1)
	if(!is.null(x$boots)){
		lines(stk$x, stk$l+1, lty = 2)
		lines(stk$x, stk$u+1, lty = 2)
	}
}


#' @rdname plot
#' @method plot vgam.crc
#' @export
plot.vgam.crc = function(x, cont.var = "x.con.1", 
		selection = NULL,
                main = "Set the main argument for the title",
                xlab = "x",
                label.offset = 0.01,
                text.size = 0.7,
                x.range = NULL, 
		subtitle = NULL,
                padj.adj = -0.3,
		lr = 1,
		lr.global = NULL,
		ylim = NULL, ...
	)
{
	xax = x$dat[,cont.var]
 	y = x$hpi
	# order y by average capture pattern frequency
	ynam = colnames(y); y = data.frame(y); names(y) = ynam
	y = y[, rev(order(colSums(y)))]
	Y = names(y)
 	ny = length(Y)
	stk = stackydens(y)
	stk$x = xax
	stk$pi0 = x$dat$pi0
	top = 1+max(stk$pi0)
	mct = x$dat$mct
	if(!is.null(x$boots)){
		stk[, c("l", "u")] = t(apply(x$boots$loc.est, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))/mct
		top = 1+max(stk$u)
	}
	if(!is.null(selection)) {
		keepers = which(rowSums(x$dat[,selection, drop = FALSE]) == length(selection))
		stk = stk[keepers,]
	}
	stk = stk[order(stk$x),]
	nx = nrow(stk)
	if(is.null(ylim)) ylim = c(0, top + 0.1)
	if(is.null(x.range)) x.range = diff(range(stk$x))
	# initialize plot
	plot(stk$x, stk[,1], , bty = "n", ylim = ylim,
		xlim = c(min(stk$x)-0.1*x.range, max(stk$x)+0.1*x.range),
		type = "n", col = 1, ylab = "Stacked relative frequency", 
		xlab = xlab, main = main, ...)
	# put in the curves
	for(k in 1:ny) lines(stk$x, stk[,k])
	# put labels on each curve
	for(k in ny:1){
		lr = -1*lr; if(!is.null(lr.global)) lr = lr.global
		if(lr > 0){
		this.x = max(stk$x) + label.offset*(x.range)
		text(Y[k], x = this.x, y = stk[nx,k], cex = text.size, adj = 0)
		}else{
		this.x= min(stk$x) - label.offset*(x.range)
		text(Y[k], x = this.x, y = stk[1,k], cex = text.size, adj = 1)
		}
	}#end for
	# add imputation
	lines(stk$x, stk$pi0+1, lwd = 2) # put in the final conditional curve
	k = nchar(names(stk)[1])
	text(paste(rep("0",k), collapse = ""), x = min(stk$x) - label.offset*(x.range), 
		y = 1 + stk$pi0[1], cex = text.size, adj = 1)
	if(!is.null(subtitle)) mtext(subtitle, padj=padj.adj, cex=text.size*1.1)
	if(!is.null(x$boots)){
		lines(stk$x, stk$l+1, lty = 2)
		lines(stk$x, stk$u+1, lty = 2)
	}
}
