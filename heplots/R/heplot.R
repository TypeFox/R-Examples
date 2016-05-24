# Change log
# last modified 23 January 2007 by J. Fox
# last modified 14 May 2007 by M. Friendly -- return xlim, ylim
# last modified 18 May 2007 by M. Friendly -- fix xlim, ylim return when !add
# last modified 20 May 2007 by M. Friendly -- pass ... to text
# last modified 23 May 2007 by J. Fox -- add ... to call to points()
# last modified 22 Oct 2007 by M. Friendly
# -- moved lambda.crit to utility.R
# -- added he.rep to handle common task of repeating HE argument values
#  13 Apr 2009 by M. Friendly -- fix label.ellipse
#  15 Apr 2009 by M. Friendly -- added axes= to fix warnings from pairs.mlm
#  24 Dec 2009 by M. Friendly -- added idate=, idesign=, icontrasts, iterm for repeated measures
#  26 Dec 2009 by M. Friendly -- workaround for car::Anova buglet
#  27 Dec 2009 by M. Friendly -- made it work for designs with no between effects
#  28 Dec 2009 by M. Friendly -- made it work with car 2.0 for doubly multivariate
#  10 Jan 2010 by M. Friendly -- merged with heplot.mlm.R
#  23 Jul 2010 by M. Friendly -- return radius
#  05 Nov 2010 by M. Friendly 
# -- added fill= and fill.alpha for filled ellipses
# -- replaced lines() with polygon() for H and E ellipses
# -- calculate H.rank to distinguish degenerate ellipses
# -- added last() to utility.R
# -- added err.label to allow changing label for Error ellipse
# -- changed default colors from palette()[-1] to a better collection, also allowing options("heplot.colors")
#  15 Jan 2013 by M. Friendly
# -- replaced internal label.ellipse with separate function; added label.pos= argument
#  22 Feb 2013
# -- added ... to label.ellipse to be able to pass cex=

`heplot` <-
		function(mod, ...) UseMethod("heplot")

`heplot.mlm` <-
		function ( 
				mod,           # an mlm object
				terms,         # vector of terms to plot H ellipses
				hypotheses,    # list of linear hypotheses for which to plot H ellipses
				term.labels=TRUE,  # TRUE, FALSE or a vector of term labels of length(terms)
				hyp.labels=TRUE,   # as above for term.labels
				err.label="Error",
				label.pos=NULL,    # label positions: NULL or 0:4
				variables=1:2,     # x,y variables for the plot [variable names or numbers]
				error.ellipse=!add,
				factor.means=!add,
				grand.mean=!add,
				remove.intercept=TRUE,
				type=c("II", "III", "2", "3"),
				idata=NULL,
				idesign=NULL,
				icontrasts=c("contr.sum", "contr.poly"),
				imatrix=NULL,
				iterm=NULL,
				markH0=!is.null(iterm),
				manova,        # an optional Anova.mlm object
				size=c("evidence", "effect.size"),
				level=0.68,
				alpha=0.05,
				segments=40,   # line segments in each ellipse
				center.pch="+",   # doesn't have to be an argument
				center.cex=2,
				col=getOption("heplot.colors", c("red", "blue", "black", "darkgreen", "darkcyan","magenta", "brown","darkgray")),
				# colors for H matrices, E matrix
				lty=2:1,
				lwd=1:2,
				fill=FALSE,         ## whether to draw filled ellipses (vectorized)
				fill.alpha=0.3,     ## alpha transparency for filled ellipses
				xlab,
				ylab,
				main="",
				xlim,           # min/max for X (override internal min/max calc) 
				ylim,
				axes=TRUE,      # whether to draw the axes
				offset.axes,    # if specified, the proportion by which to expand the axes on each end (e.g., .05)
				add=FALSE,      # add to existing plot?
				verbose=FALSE,
				warn.rank=FALSE,  
				...) {
	ell <- function(center, shape, radius) {
		angles <- (0:segments)*2*pi/segments
		circle <- radius * cbind( cos(angles), sin(angles))
		if (!warn.rank){
			warn <- options(warn=-1)
			on.exit(options(warn))
		}
		Q <- chol(shape, pivot=TRUE)
		order <- order(attr(Q, "pivot"))
		t( c(center) + t( circle %*% Q[,order]))
	}
#	label.ellipse <- function(ellipse, label, col){
#		if (cor(ellipse)[1,2] >= 0){
#			index <- which.max(ellipse[,2])
#			x <- ellipse[index, 1] + 0.5 * strwidth(label)  # was: "A"
#			y <- ellipse[index, 2] + 0.5 *strheight("A")
#			adj <- c(1, 0) 
#		}
#		else {
#			index <- which.min(ellipse[,2])
#			x <- ellipse[index, 1] - 0.5 * strwidth(label)  # was: "A"
#			y <- ellipse[index, 2] - 0.5 * strheight("A")
#			adj <- c(0, 1) 
#		}
#		text(x, y, label, adj=adj, xpd=TRUE, col=col, ...)
#	}
#	last <- function(x) {x[length(x)]}
	
	#if (!require(car)) stop("car package is required.")
	# avoid deprecated warnings from car
	if (packageDescription("car")[["Version"]] >= 2) linear.hypothesis <- linearHypothesis
	type <- match.arg(type)
	size <- match.arg(size)
	data <- model.frame(mod)
	
	if (missing(manova)) {
		if (is.null(imatrix)) {
			manova <- Anova(mod, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts)
		}
		else {
			if (packageDescription("car")[["Version"]] >= 2)
				manova <- Anova(mod, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)
			else stop("imatrix argument requires car 2.0-0 or later")
		} 
	}   
	if (verbose) print(manova)
	
	if (is.null(idata) && is.null(imatrix)) {
		Y <- model.response(data) 
		SSPE <- manova$SSPE
	} 
	else {
		if (is.null(iterm)) stop("Must specify a within-S iterm for repeated measures designs" )
		### DONE::car -- workaround for car::Anova.mlm bug: no names assigned to $P component
		if (is.null(names(manova$P))) names(manova$P) <- names(manova$SSPE)
		Y <- model.response(data) %*% manova$P[[iterm]]
		SSPE <- manova$SSPE[[iterm]]
	}   
	
	if (!is.null(rownames(SSPE))) {response.names <- rownames(SSPE)}
	else {response.names <- paste("V.", 1:nrow(SSPE), sep="")}
	p <- length(response.names)
	
	if (!is.numeric(variables)) {
		vars <- variables
		variables <- match(vars, response.names)
		check <- is.na(variables)
		if (any(check)) stop(paste(vars[check], collapse=", "), 
					" not among response variables.") 
	}
	else {
		if (any (variables > length(response.names))) stop("There are only ", 
					length(response.names), " response variables.")
		vars <- response.names[variables]
	}
	if (length(variables) != 2) {
		extra <- if (length(variables) == 1) 'heplot1d()' else 
				if (length(variables) == 3) 'heplot3d()' else 'pairs()'
		stop(paste("You may only plot 2 response variables. Use", extra))
	}
	
	if (missing(terms) || (is.logical(terms) && terms)) {
		terms <- manova$terms
# FIXME:  This does mot work if the between-S design includes only an intercept 
# FIXME: && terms="(Intercept)" is specified
		if (!is.null(iterm)) {
#			if (terms=="(Intercept)")  terms <- iterm else 
			terms <- terms[grep(iterm, terms)]   ## only include those involving iterm
		}
		if (remove.intercept) terms <- terms[terms != "(Intercept)"]
	}
	n.terms <- if (!is.logical(terms)) length(terms) else 0 
	# note: if logical here, necessarily FALSE
	n.hyp <- if (missing(hypotheses)) 0 else length(hypotheses)
	n.ell <- n.terms + n.hyp
	if (n.ell == 0) stop("Nothing to plot.")
	
	Y <- Y[,vars] 
	gmean <- if (missing(data))  c(0,0) 
			else colMeans(Y)
	if (missing(xlab)) xlab <- vars[1]
	if (missing(ylab)) ylab <- vars[2] 
	dfe <- manova$error.df
	scale <- 1/dfe
	radius <- sqrt(2 * qf(level, 2, dfe))
	
	# assign colors and line styles
	col <- he.rep(col, n.ell) 
	lty <- he.rep(lty, n.ell)
	lwd <- he.rep(lwd, n.ell)
	# handle filled ellipses
	fill <- he.rep(fill, n.ell)
	fill.alpha <- he.rep(fill.alpha, n.ell)
	fill.col <- trans.colors(col, fill.alpha)
	label.pos <- he.rep(label.pos, n.ell)
	# TODO:  take account of rank=1?
	fill.col <- ifelse(fill, fill.col, NA)
	E.col<- last(col)
	
	H.ellipse <- as.list(rep(0, n.ell))
	# keep track of ranks to distinguish degenerate ellipses
	H.rank <- rep(0, n.ell)
	
	if (n.terms > 0) for (term in 1:n.terms){
			term.name <- terms[term]
			H <- manova$SSP[[term.name]]
			if (!(all(variables %in% 1:nrow(H)))) {
				warning(paste("Skipping H term ", term.name, "(size: ", nrow(H), ")", sep=""))
				next
			}
			H <- H[variables, variables]
			dfh <- manova$df[term.name]
			factor <- if (size == "evidence") lambda.crit(alpha, p, dfh, dfe) else 1
			H <- H * scale/factor
			if (verbose){
				cat(term.name, " H matrix (", dfh, " df):\n")
				print(H)
			}
			H.ellipse[[term]] <- ell(gmean, H, radius)
			H.rank[term] <- qr(H)$rank
		}
	if (n.hyp > 0) for (hyp in 1:n.hyp){
			lh <- linear.hypothesis(mod, hypotheses[[hyp]])
			H <- lh$SSPH[variables, variables]
			dfh <- lh$df
			factor <- if (size == "evidence") lambda.crit(alpha, p, dfh, dfe) else 1
			H <- H * scale/factor
			if (verbose){
				cat("\n\n Linear hypothesis: ", names(hypotheses)[[hyp]], "\n") 
				print(lh)
			}
			H.ellipse[[n.terms + hyp]] <- ell(gmean, H, radius)
		}
	E <- SSPE
	E <- E[variables, variables]
	E <- E * scale[1]
	E.ellipse <- ell(gmean, E, radius)
	H.ellipse$E <- E.ellipse     
	if (!add){
		max <- apply(sapply(H.ellipse, function(X) apply(X, 2, max)), 1, max)
		min <- apply(sapply(H.ellipse, function(X) apply(X, 2, min)), 1, min)
		factors <- data[, sapply(data, is.factor), drop=FALSE]
		if (!is.logical(factor.means)){
			factor.names <- colnames(factors) 
			which <- match(factor.means, factor.names)
			check <- is.na(which)
			if (any(check)) stop(paste(factor.means[check], collapse=", "), 
						" not among factors.")
			factors <- factors[, which, drop=FALSE]
		}
		if (!is.logical(factor.means) || factor.means){   
			for (fac in factors){
				means <- aggregate(Y, list(fac), mean)
				min[1] <- min(min[1], means[,2])
				max[1] <- max(max[1], means[,2])
				min[2] <- min(min[2], means[,3])
				max[2] <- max(max[2], means[,3])
			}
		}
		if (!missing(offset.axes)){
			range <- max - min
			min <- min - offset.axes*range
			max <- max + offset.axes*range
		}
		xlim <- if(missing(xlim)) c(min[1], max[1]) else xlim
		ylim <- if(missing(ylim)) c(min[2], max[2]) else ylim
		plot(xlim, ylim,  type = "n", xlab=xlab, ylab=ylab, main=main, axes=axes, ...)
	}
	# no longer need H.ellipse$E, since we return it separately
	H.ellipse$E <- NULL
	if (grand.mean) 
		points(gmean[1], gmean[2], pch=center.pch, cex=center.cex, col="black", xpd=TRUE)
	if (error.ellipse){
#		lines(E.ellipse, col=E.col, lty=lty[length(lty)], lwd=lwd[length(lwd)])
		polygon(E.ellipse, col=last(fill.col), border=last(col), lty=last(lty), lwd=last(lwd))
		label.ellipse(E.ellipse, err.label, col=last(col), label.pos=last(label.pos), ...)
	}
	term.labels <- if (n.terms == 0) NULL
			else if (!is.logical(term.labels)) term.labels
			else if (term.labels) terms else rep("", n.terms)  
	if (n.terms > 0) for (term in 1:n.terms){
#			lines(H.ellipse[[term]], col=col[term], lty=lty[term], lwd=lwd[term])
			# TODO: avoid polygon if rank=1 ???
			polygon(H.ellipse[[term]], col=fill.col[term], border=col[term],  lty=lty[term], lwd=lwd[term])
			label.ellipse(H.ellipse[[term]], term.labels[term], col=col[term], label.pos=label.pos[term], ...) 
		}   
	hyp.labels <- if (n.hyp == 0) NULL
			else if (!is.logical(hyp.labels)) hyp.labels
			else if (hyp.labels) names(hypotheses) else rep("", n.hyp)  
	if (n.hyp > 0) for (hyp in 1:n.hyp){
			ell <- n.terms + hyp
#			lines(H.ellipse[[ell]], col=col[ell], lty=lty[ell], lwd=lwd[ell])
			polygon(H.ellipse[[ell]], col=fill.col[ell], border=col[ell],  lty=lty[ell], lwd=lwd[ell])
			label.ellipse(H.ellipse[[ell]], hyp.labels[hyp], col=col[ell], label.pos=label.pos[ell], ...)
		}
	if (!add && (!is.logical(factor.means) || factor.means)){
		for (fac in factors){
			means <- aggregate(Y, list(fac), mean)
			points(means[,2], means[,3], pch=16, xpd=TRUE, ...)
			text(means[,2], means[,3], labels=as.character(means[,1]), pos=3, xpd=TRUE, ...)
		}
	}
	
	if(is.logical(markH0) && markH0) mark.H0()
	else if (is.list(markH0)) do.call(mark.H0, markH0)
	
	names(H.ellipse) <- c(if (n.terms > 0) term.labels, if (n.hyp > 0) hyp.labels)
	result <- if (!add) list(H=H.ellipse, E=E.ellipse, center=gmean, xlim=xlim, ylim=ylim, radius=radius)
			else list(H=H.ellipse, E=E.ellipse, center=gmean, radius=radius)
	class(result) <- "heplot"
	invisible(result)
}
