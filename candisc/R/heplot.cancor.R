# HE plot for a cancor object

heplot.cancor <- function (
	mod,		         # output object from cancor
	which=1:2,       # canonical dimensions to plot
	scale,           # scale factor(s) for variable vectors in can space
	asp=NA,           # aspect ratio, to ensure equal units?
	var.vectors = "Y", # which variable vectors to show?
	var.col=c("blue", "darkgreen"),  # colors for Y and X variable vectors and labels
	var.lwd=par("lwd"),
	var.cex=par("cex"),
	var.xpd=TRUE,     # was: par("xpd"),
	prefix = "Ycan",  # prefix for labels of canonical dimensions
	suffix = TRUE,   # add label suffix with can % ?
	terms=TRUE,  # terms to be plotted in canonical space / TRUE=all
	...              # extra args passed to heplot
	) {

  if (!inherits(mod, "cancor")) stop("Not a cancor object")
	if (mod$ndim < 2 || length(which)==1) {
		# using stop() here would terminate heplot.candiscList
	   message("Can't do a 1 dimensional canonical HE plot")
#	   plot(mod, which=which, var.col=var.col, var.lwd=var.lwd, prefix=prefix, suffix=suffix, ...) 
	   return()
	}

	Yvars <- mod$names$Y
	scores <- data.frame(mod$scores$X, mod$scores$Y)
	scores <- data.frame(scores, mod$X)   # append X variables
	Xcoef <- mod$coef$X
	Ycoef <- mod$coef$Y
	Ycan <- colnames(Ycoef)

	canr <- mod$cancor
  lambda <- canr^2 / (1-canr^2)
  pct = 100*lambda / sum(lambda)

  if ((is.logical(terms) && terms) || terms=="X") {
  	terms <- mod$names$X
	}
	# allow plotting the Xcan variables
	else if (length(terms)==1 && terms=="Xcan") terms=colnames(Xcoef)

	# make sure that all terms are available	
	if (!all(terms %in% colnames(scores))) {
			stop(paste(setdiff(terms, colnames(scores) ), "are not among the available variables"))
		}

##   Construct the model formula to fit mod$Yscores ~ Xscores in original lm()
##   using the mod$scores data.frame
#browser()
  txt <- paste( "lm( cbind(",
              paste(Ycan, collapse = ","),
              ") ~ ",
              paste( terms, collapse = "+"), ", data=scores)" )
  can.mod <- eval(parse(text=txt))
  
##   Construct labels for canonical variables
	canvar <- Ycan[which]   # names of canonical variables to plot
	if (is.logical(suffix) & suffix)
		suffix <- paste( " (", round(pct[which],1), "%)", sep="" ) else suffix <- NULL
	canlab <- paste(prefix, which, suffix, sep="")

  ellipses <- heplot(can.mod, terms=terms, 
  		xlab=canlab[1], ylab=canlab[2],  ...)
  
	struc <- mod$structure
  Xstructure <- struc$X.yscores[,which]
  Ystructure <- struc$Y.yscores[,which]
  vec <- rbind(
  	if ("Y" %in% var.vectors) Ystructure else NULL,
  	if ("X" %in% var.vectors) Xstructure else NULL)

	maxrms <- function(x) { max(sqrt(apply(x^2, 1, sum))) }
	if (missing(scale)) {
		vecrange <- range(vec)
		ellrange <- lapply(ellipses, range)
		vecmax <- maxrms(vec)
		ellmax <- max( maxrms(ellipses$E), unlist(lapply(ellipses$H, maxrms)) )
		scale <- floor(  0.9 * ellmax / vecmax )
		cat("Vector scale factor set to ", scale, "\n")
	}
  
  if ("Y" %in% var.vectors) vectors(Ystructure, scale=scale, col=var.col[1], lwd=var.lwd, cex=var.cex, xpd=var.xpd)
  if ("X" %in% var.vectors) vectors(Xstructure, scale=scale, col=var.col[2], lwd=var.lwd, cex=var.cex, xpd=var.xpd)

	invisible(ellipses)
	

}
