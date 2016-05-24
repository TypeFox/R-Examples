## 1D representation of an HE plot
# Initial version 17-Apr-2009
# Fixed buglet with hyp.labels 8-Dec-2009
# last modified 1 Jan 2010 by M. Friendly -- added idate=, idesign=, icontrasts, iterm for repeated measures


`heplot1d` <-
function(mod, ...) UseMethod("heplot1d")


`heplot1d.mlm` <-
		function ( 
				mod,           # an mlm object
				terms,         # vector of terms to plot H ellipses
				hypotheses,    # list of linear hypotheses for which to plot H ellipses
				term.labels=TRUE,  # TRUE, FALSE or a vector of term labels of length(terms)
				hyp.labels=TRUE,   # as above for term.labels
				variables=1,       # x,y variables for the plot [variable names or numbers]
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
				manova,        # an optional Anova.mlm object
				size=c("evidence", "effect.size"),
				level=0.68,
				alpha=0.05,
				center.pch="|",   # doesn't have to be an argument
				col=getOption("heplot.colors", c("red", "blue", "black", "darkgreen", "darkcyan","magenta", "brown","darkgray")),
				# colors for H matrices, E matrix
				lty=2:1,
				lwd=1:2,
				xlab,
				main="",
				xlim,           # min/max for X (override internal min/max calc) 
				axes=TRUE,      # whether to draw the axes
				offset.axes=0.1,    # proportion by which to expand the axes on each end (e.g., .05)
				add=FALSE,      # add to existing plot?
				verbose=FALSE,
				...) {

	ell1d <- function(center, shape, radius) {
		circle <- radius * c(-0.5, 0.5)
		center + sqrt(shape) * circle 
	}
	F.crit <- function(alpha, p, dfh, dfe) {
		(dfh/dfe) * qf(alpha, dfh, dfe, lower.tail=FALSE)
	}

	#if (!require(car)) stop("car package is required.")
	if (car2 <- packageDescription("car")[["Version"]] >= 2) linear.hypothesis <- linearHypothesis
	type <- match.arg(type)
	size <- match.arg(size)
	data <- model.frame(mod)

#	if (missing(manova)) manova <- Anova(mod, type=type)    
	if (missing(manova)) {
		if (is.null(imatrix)) {
			manova <- Anova(mod, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts)
		}
		else {
			if (car2)
				manova <- Anova(mod, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)
			else stop("imatrix argument requires car 2.0-0 or later")
		} 
	}   
	
#	if (verbose) print(manova)    

	if (is.null(idata) && is.null(imatrix)) {
		Y <- model.response(data) 
		SSPE <- manova$SSPE
	} 
	else {
		if (is.null(iterm)) stop("Must specify a within-S iterm for repeated measures designs" )
		### FIXME::car -- workaround for car::Anova.mlm bug: no names assigned to $P component
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

### Allow for more than one variable?
	if (length(variables) != 1) {
#		stop(paste("You may only plot 1 response variable."))
		extra <- if (length(variables) == 2) 'heplot()' else 
			if (length(variables) == 3) 'heplot3d()' else 'pairs()'
		stop(paste("You may only plot 1 response variable. Try", extra))
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
	gmean <- if (missing(data))  0 else mean(Y)
#			else colMeans(Y)

	if (missing(xlab)) xlab <- vars[1]
	dfe <- manova$error.df
	scale <- 1/dfe
#	radius <- sqrt(2 * qf(level, 2, dfe))
	radius <- sqrt(    qf(level, 1, dfe))

	# assign colors and line styles
	col <- he.rep(col, n.ell); E.col<- col[length(col)]
	lty <- he.rep(lty, n.ell)
	lwd <- he.rep(lwd, n.ell)

	# plot the 1D representations of the terms on equally spaced lines
	yvals <- 1:n.ell

	H.ellipse <- as.list(rep(0, n.ell))
	if (n.terms > 0) for (term in 1:n.terms){
			term.name <- terms[term]
			H <- manova$SSP[[term.name]]
			H <- H[variables, variables]
			dfh <- manova$df[term.name]
			factor <- if (size == "evidence") F.crit(alpha, p, dfh, dfe) else 1
			H <- H * scale/factor
			if (verbose){
				cat(term.name, " H matrix (", dfh, " df):\n")
				print(H)
			}
			H.ellipse[[term]] <- ell1d(gmean, H, radius)
			if(verbose) {cat(term.name, "H range:\n"); print(H.ellipse[[term]])}
		}
	if (n.hyp > 0) for (hyp in 1:n.hyp){
#			lh <- linear.hypothesis(mod, hypotheses[[hyp]])
			lh <- linearHypothesis(mod, hypotheses[[hyp]])
			H <- lh$SSPH[variables, variables]
			dfh <- lh$df
			factor <- if (size == "evidence") F.crit(alpha, p, dfh, dfe) else 1
			H <- H * scale/factor
			if (verbose){
				cat("\n\n Linear hypothesis: ", names(hypotheses)[[hyp]], "\n") 
				print(lh)
			}
			H.ellipse[[n.terms + hyp]] <- ell1d(gmean, H, radius)
		}
	E <- SSPE
	E <- E[variables, variables]
	E <- E * scale[1]
	E.ellipse <- ell1d(gmean, E, radius)
	H.ellipse$E <- E.ellipse     

	if (!add){
#		max <- apply(sapply(H.ellipse, function(X) apply(X, 2, max)), 1, max)
#		min <- apply(sapply(H.ellipse, function(X) apply(X, 2, min)), 1, min)
		max <- max(sapply(H.ellipse, max))
		min <- min(sapply(H.ellipse, min))

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
				min <- min(min, means[,2])
				max <- max(max, means[,2])
			}
		}
		if (!missing(offset.axes)){
			range <- max - min
			min <- min - offset.axes*range
			max <- max + offset.axes*range
		}
		xlim <- if(missing(xlim)) c(min[1], max[1]) else xlim
#		ylim <- if(missing(ylim)) c(min[2], max[2]) else ylim
		plot(xlim, range(yvals),  type = "n", xlab=xlab, 
				ylab=if (n.hyp>0) "Terms and Hypotheses" else "Terms", yaxt="n",
				main=main, axes=axes, ...)
	}

	H.ellipse$E <- NULL
	if (error.ellipse){
#		lines(E.ellipse, col=E.col, lty=lty[length(lty)], lwd=lwd[length(lwd)])
		rect(E.ellipse[1], 0.5, E.ellipse[2], n.ell+0.5, 
		col=E.col, lty=lty[length(lty)], lwd=lwd[length(lwd)], border=NA)
#		label.ellipse(E.ellipse, "Error", col=E.col)
	}
	if (grand.mean) 
		points(rep(gmean,n.ell), 1:n.ell, pch=center.pch, cex=2, col="black", xpd=TRUE)

	term.labels <- if (n.terms == 0) NULL
			else if (!is.logical(term.labels)) term.labels
			else if (term.labels) terms else rep("", n.terms)  
	if (n.terms > 0) for (term in 1:n.terms){
			lines(x=H.ellipse[[term]], y=rep(term,2), col=col[term], lty=lty[term], lwd=lwd[term])
#			label.ellipse(H.ellipse[[term]], term.labels[term], col=col[term]) 
			text(xlim[1],term, term.labels[term], col=col[term], adj=c(0,0))
			term.name <- terms[term]
			means <- termMeans(mod, term.name, label.factors=FALSE)
			points(means[,vars], rep(term,nrow(means)), pch=16, xpd=TRUE, ...)
			widths <- strwidth(rownames(means))
# TODO: determin pos based on whether there is overlap of labels
			pos <- rep(c(1,3),length=nrow(means))
			text(means[,vars], rep(term,nrow(means)), labels=rownames(means), 
					pos=pos, xpd=TRUE, ...)
			if (verbose){
				cat("\n",term.name, " means:\n")
				print(means[,vars,drop=FALSE])
			}
		}   
	hyp.labels <- if (n.hyp == 0) NULL
			else if (!is.logical(hyp.labels)) hyp.labels
			else if (hyp.labels) names(hypotheses) else rep("", n.hyp)  
	if (n.hyp > 0) for (hyp in 1:n.hyp){
			term <- n.terms + hyp
			lines(x=H.ellipse[[term]], y=rep(term,2), col=col[term], lty=lty[term], lwd=lwd[term])
#			label.ellipse(H.ellipse[[term]], hyp.labels[hyp], col=col[term])
			text(xlim[1],term, hyp.labels[term], col=col[term], adj=c(0,0))
		}

#	if (!add && (!is.logical(factor.means) || factor.means)){
#		line <- 0
#		for (fac in factors){
#			line <- line+1
#			means <- aggregate(Y, list(fac), mean)
#			if (verbose){
#				cat(colnames(factors)[fac], " means:\n")
#				print(means)
#			}
#			points(means[,2], rep(line,nrow(means)), pch=16, xpd=TRUE, ...)
#			text(means[,2], rep(line,nrow(means)), labels=as.character(means[,1]), 
#					pos=rep(c(1,3),length=nrow(means)), xpd=TRUE, ...)
#		}
#	}

	names(H.ellipse) <- c(if (n.terms > 0) term.labels, if (n.hyp > 0) hyp.labels)
	result <- if (!add) list(H=H.ellipse, E=E.ellipse, xlim=xlim)	else list(H=H.ellipse, E=E.ellipse)
	class(result) <- "heplot1d"
	invisible(result)

}

