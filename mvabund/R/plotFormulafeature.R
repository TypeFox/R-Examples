
################################################################################
# plot formulas, should not be used directly, use plot.mvformula instead  #
# similar to plot.formula, but has some additional features:                   #
# - shift of overlapping points possible                                       #
# - the axes can be specified with vectors yaxis.ticks and yaxis.labs and      #
# lists xaxis.ticks, xaxis.labs                                                #
################################################################################

plotFormulafeature <- function( formula, 
				data = parent.frame(),
				subset, 
				ylab = varnames[response],
				xlab=NULL, 
				ask = TRUE, 
				las=1, 
				yaxis.ticks=NULL, 
				yaxis.labs=NULL,
				xaxis.ticks=NULL, 
				xaxis.labs=NULL, 
				fg="grey",
				cex.xaxis=0.9, 
				cex.yaxis=0.6, 
				overall.main="", 
				main="",
				axes=!any(!is.null(c(yaxis.ticks, yaxis.labs, xaxis.ticks, xaxis.labs))) , 
				xvar.subset=NULL, 
				all.labels=FALSE,
				 ... ) {

mfrows <- par("mfrow")
colmains <- par("col.main")

pw <- mfrows[1]*mfrows[2]
if (pw ==1 & overall.main!="") main <- paste("\n \n",main)
	
m <- match.call(expand.dots = FALSE)
	
if (is.matrix(eval(m$data, parent.frame()))) 

m$data <- as.data.frame(data)
dots <- m$...
	
dots <- lapply(dots, eval, data, parent.frame())
	
if ("type" %in% names(dots)) type <- dots$type
else type <- "p"
	
if ("col" %in% names(dots)) {
	dcol <- dots$col
	dots$col <- NULL

#	if(type=="bx") dcolbx <- dcol
#	else dcolbx <- NULL

} else {
	dcol <- "black"
#	dcolbx <- NULL
}

if ("border" %in% names(dots)) {
	bordcol <- dots$border
	dots$border <- NULL
} else	bordcol <- par("fg")
	
pch <- dots$pch

m$ylab <-  m$xlab <- m$... <- m$ask <- m$las <- m$yaxis.ticks <- m$yaxis.labs <- m$xaxis.ticks <- m$xaxis.labs <- m$fg <- m$cex.xaxis <- m$cex.yaxis  <- m$axes <- m$main <- m$overall.main <- m$xvar.subset <-m$all.labels <- NULL

subset.expr <- m$subset

m$subset <- NULL
m[[1]] <- as.name("model.frame")
m <- as.call(c(as.list(m), list(na.action = NULL)))   
mf <- eval(m, parent.frame())
model.mat <- model.matrix(formula, data= data)

is.interact <- extend.x.formula(formula, return.interaction=TRUE,extend.term=TRUE)$is.interaction

# Exclude Interactions with factors, as it would involve too much hassle to get x.
if(any(is.interact)){
	terms.foo <- attr(model.frame(formula, data=data), "terms")

	facs <- attr(terms.foo,"factors")[,is.interact, drop = FALSE]
	datClass <- attr(terms.foo,"dataClasses")

	for(faci in 1:NCOL(facs)){
		if(any(datClass[facs[,faci]>0] == "factor"))
			stop("Formulas with interactions that include factors cannot yet be plotted")
	}
}

if (!missing(subset)) {
	s <- eval(subset.expr, data, parent.frame())
	l <- nrow(mf)
	dosub <- function(x) if (length(x) == l) x[s]
			     else x
	dots <- lapply(dots, dosub)
	mf <- mf[s, , drop=FALSE ]
	model.mat <- model.mat[s,, drop=FALSE]
}
	
horizontal <- FALSE
if ("horizontal" %in% names(dots)) 
	horizontal <- dots[["horizontal"]]
response <- attr(attr(mf, "terms"), "response")

if (response) {
	varnames <- names(mf)
	y <- mf[[response]]
	funname <- NULL
	
	if (is.object(y)) {
		found <- FALSE
		if (is.mvabund(y))
			y <- unabund(y)
	
		for (j in class(y)) {
			funname <- paste("plot.", j, sep = "")
			if (exists(funname)) {
				found <- TRUE
				break
			}
		}

		if (!found ) funname <- NULL
	}

	if (is.null(funname)) {funname <- "plot"}	
	if (length(varnames) > 2) {
		opar <- par(ask = ask)
		on.exit(par(opar))
	}

	xn <- attr(terms(formula), "term.labels") # varnames[-response]
	if (length(xn) > 0) {
		if (! is.null(xlab)){ 
			if(length(xlab)==1) {
				xlab<-rep(xlab,times =length(xn))
			} else if (length(xlab)!=length(xn)) {
				stop("'xlab' does not have the appropriate lenght")
			}
		}	  
		
		j <- 0
		jj <- 1:length(xn)
		if (!is.null(xvar.subset)) {
			if(is.logical(xvar.subset) & any(!is.na(xvar.subset)))
			
			var.subset <- which(xvar.subset[!is.na(xvar.subset)])
			
			if(min(xvar.subset)<1 | max(xvar.subset)>length(xn))
			stop("'xvar.subset' is invalid, 'xvar.subset' must be a vector of integers between 1 and ",length(xn))
			
			xn <- xn[xvar.subset]
			jj <- jj[xvar.subset]
		}

		for (i in xn) {
			j <- j+1
			xl <- if (is.null(xlab)) i
			      else xlab[jj[j]]

			if (all.labels | i==xn[1] ) yl <- ylab
			else yl <- ""	

			# For Interactions get x.data via the model.matrix, as model.frame does
			# represent them differently.
			if(is.interact[jj[j]]) { 
				x.data <-  model.mat[ , colnames(model.mat) == i]
			} else  x.data <- mf[[i]]
				
#			if (type=="bx") x.data <- factor(x.data)
	
			if (horizontal && is.factor(x.data)) {
				yl <- xl
				xl <- ylab
			}

			if(is.list(dcol)) dcoli <- dcol[[jj[j]]] 
			else dcoli <- dcol

#			if(is.list(dcolbx)) dcolbxi <- dcolbx[[jj[j]]] 
#			else dcolbxi <- dcolbx

			if(is.list(pch)) dots$pch <- pch[[jj[j] ]]
			


			if (is.factor(x.data)) {
cat("plot 1\n")			
				do.call(funname, c(list(x.data, y, ylab = "", xlab = "", col=dcol,
						border=bordcol) , main=main, axes =axes,las=las, fg=fg, dots))
				title(xlab=xl, ylab=yl, cex.lab=1.2)
			} else { 
cat("plot 2\n")
				do.call(funname, c(list(x.data, y, ylab = "", xlab = "", col=dcoli),
						main=main, axes =axes,las=las, fg=fg, dots))
				title(xlab=xl, ylab=yl, cex.lab=1.2)
			}

			if (!missing(axes) ) {
				if(!axes) {
					box(col=fg)
					axis( side=2, at=yaxis.ticks, labels=yaxis.labs, las=las,
						col=fg, cex.axis=cex.yaxis)
					axis( side=1, at=xaxis.ticks[[i]], labels=xaxis.labs[[i]], las=las,
						col=fg, cex.axis=cex.xaxis)
				}
			}

		}
        	mtext(overall.main, outer = TRUE, cex = par("cex.main"), col=colmains, line=1)

	} else {
		if (type=="bx") {
#			do.call("boxplot", c(list(y, ylab = ylab, col=dcolbx, border=bordcol),
#					main=main, axes =axes,las=las, fg=fg, dots))
		} else {
cat("plot 3\n")
			do.call(funname, c(list(y, ylab = ylab, col=dcol), main=main,
					axes =axes, las=las, fg=fg, dots))
		}
		
		if (!missing(axes) ) {
			if(!axes)  {
				box(col=fg)
				axis( side=2, at=yaxis.ticks, labels=yaxis.labs, las=las,
					col=fg, cex.axis=cex.yaxis)
				axis( side=1, las=las, col=fg, cex.axis=cex.xaxis)
			} 
		}
		
		mtext(overall.main, outer = TRUE, cex = par("cex.main"), col=colmains) 
				
	}

} else {
	stop("'formula' has no response. Use standard plotting functions") 
}

}

