# HE plot for a cancor object

heplot3d.cancor <- function (
	mod,		         # output object from cancor
	which=1:3,       # canonical dimensions to plot
	scale,           # scale factor(s) for variable vectors in can space
	asp="iso",           # aspect ratio, to ensure equal units
	var.vectors = "Y", # which variable vectors to show? [not yet implemented]
	var.col=c("blue", "darkgreen"),  # colors for Y and X variable vectors and labels
	var.lwd=par("lwd"),
	var.cex=par("cex"),
	var.xpd=NA,       # not used
	prefix = "Ycan",  # prefix for labels of canonical dimensions
	suffix = FALSE,   # add label suffix with can % ?
	terms=TRUE,  # terms to be plotted in canonical space / TRUE=all
	...              # extra args passed to heplot
	) {

  if (!inherits(mod, "cancor")) stop("Not a cancor object")
	if (mod$ndim < 3 || length(which) < 3) {
		# using stop() here would terminate heplot.candiscList
	   message("Can't do a 3 dimensional canonical HE plot")
#	   plot(mod, which=which, var.col=var.col, var.lwd=var.lwd, prefix=prefix, suffix=suffix, ...) 
	   return()
	}
  if (!requireNamespace("rgl")) stop("rgl is required")

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

  ellipses <- heplot3d(can.mod, terms=terms, 
  		xlab=canlab[1], ylab=canlab[2], zlab=canlab[3], ...)
  
	struc <- mod$structure
  Xstructure <- struc$X.yscores[,which]
  Ystructure <- struc$Y.yscores[,which]


# TODO: calculate scale factor(s)
  structure <- Ystructure
	maxrms <- function(x) { max(sqrt(apply(x^2, 1, sum))) }
	if (missing(scale)) {
		vecmax <- maxrms(structure)
		vecrange <- range(structure)
		ellrange <- lapply(ellipses, range)
		vecmax <- maxrms(structure)
		# get bbox of the 3d plot
        bbox <- matrix(rgl::par3d("bbox"),3,2,byrow=TRUE)
#    TODO: calculate scale so that vectors reasonably fill the plot
		scale <- 3
		cat("Vector scale factor set to ", scale, "\n")
#	  browser()
	}

# TODO: plot vectors

  cs <- scale * Ystructure
  #  can this be simplified?
  for(i in 1:nrow(cs)) {
  	rgl::lines3d( c(0, cs[i,1]),
  	         c(0, cs[i,2]),
  	         c(0, cs[i,3]), col=var.col, lwd=var.lwd)
  }
  rgl::texts3d( cs, texts=rownames(cs), col=var.col, cex=var.cex)


  if (!is.null(asp)) rgl::aspect3d(asp)
  
  invisible(ellipses)

}
