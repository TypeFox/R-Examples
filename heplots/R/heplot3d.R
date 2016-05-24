# last modified 23 January 2007 by J. Fox
# last modified 10/6/2008 8:32AM by M. Friendly
#    - added shade=, shade.alpha=, wire=
#    - fixed: grand.mean=FALSE not respected
#    - replaced sphere at grand.mean with cross3d
# last modified 10/7/2008 11:26AM by M. Friendly
#    - color means according to the color of the term
# last modified 10/15/2008 4:21PM by M. Friendly
#    - return bounding boxes of the ellipsoids
# last modified 10/28/2008 9:37AM by M. Friendly
#    - replaced rgl.texts with texts3d
#    - replaced ellipsoid() with spheres3d() for plotting means
# last modified 11/4/2008 12:50PM by M. Friendly
#    - reverted to ellipsoid() for plotting means
# last modified 11/6/2008 by M. Friendly
#    - added xlim, ylim, zlim arguments to allow expanding the plot bbox (for candisc)
# last modified 29 Dec 2009 by M. Friendly -- added idate=, idesign=, icontrasts, iterm for repeated measures
# last modified 30 Dec 2009 by M. Friendly -- debugged repeated measures
# last modified  1 Jan 2010 by M. Friendly -- debugged repeated measures again
# last modified  1 Jan 2010 by M. Friendly -- merged heplot3d.R and heplot.mlm.R
# last modified 12 Feb 2010 by M. Friendly -- fixed buglet with text3d causing rgl to crash (thx: Duncan Murdoch)
# last modified 23 Jul 2010 by M. Friendly -- return radius
# -- added err.label to allow changing label for Error ellipse
# last modified 26 Apr 2013 by M. Friendly 
# -- modified ellipsoid to reduce striation (Thx: Duncan Murdoch)
# -- changed default colors and default fill.alpha

savedvars <- new.env(parent=emptyenv())

`heplot3d` <-
function(mod, ...) UseMethod("heplot3d")

# TODO:
#  - add aspect argument (for candisc)

`heplot3d.mlm` <-
		function ( 
				mod,           # an mlm object
				terms,         # vector of terms to plot H ellipses
				hypotheses,    # list of linear hypotheses for which to plot H ellipses
				term.labels=TRUE,  # TRUE, FALSE or a list of term labels of length(terms)
				hyp.labels=TRUE,   # as above for term.labels
				err.label="Error",
				variables=1:3,     # x,y variables for the plot [variable names or numbers]
				error.ellipsoid=!add,
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
				segments=40,          # line segments in each ellipse
#				col=palette()[-1],    # colors for E matrix, H matrices
				col=getOption("heplot3d.colors", c("red", "blue", "black", "darkgreen", "darkcyan","magenta", "brown","darkgray")),
				# colors for H matrices, E matrix
				lwd=c(1, 4),          # line width for drawing ellipsoids and 1d degenerate ellipsoids
				shade=TRUE,           # use shade3d to render ellipsoids?
				shade.alpha=0.2,      # alpha transparency for shaded3d
				wire=c(TRUE,FALSE),   # use wire3d to render ellipsoids?
				bg.col=c("white", "black"),  # background colour
				fogtype=c("none", "exp2", "linear", "exp"), # fog -- for depth cueing
				fov=30,   # field of view (for perspective)
				offset=0.01, # for ellipsoid labels
				xlab,
				ylab,
				zlab,
				xlim,
				ylim,
				zlim,
				add=FALSE,      # add to existing plot?
				verbose=FALSE,
				warn.rank=FALSE,  
				...) {              
	
	ellipsoid <- function(center, shape, radius=1, label="", col, df=Inf, shade=TRUE, alpha=0.1, wire=TRUE){
		# adapted from the shapes3d demo in the rgl package and from the Rcmdr package
		# modified to return the bbox of the ellipsoid
		degvec <- seq(0, 2*pi, length=segments)
		ecoord2 <- function(p) c(cos(p[1])*sin(p[2]), sin(p[1])*sin(p[2]), cos(p[2]))
		# v <- t(apply(expand.grid(degvec,degvec), 1, ecoord2))  # modified to make smoother
		v <- t(apply(expand.grid(degvec,degvec/2), 1, ecoord2)) 
		if (!warn.rank){
			warn <- options(warn=-1)
			on.exit(options(warn))
		}
		Q <- chol(shape, pivot=TRUE)
		lwd <- if (df < 2) lwd[2] else lwd[1]
		order <- order(attr(Q, "pivot"))
		v <- center + radius * t(v %*% Q[, order])
		v <- rbind(v, rep(1,ncol(v))) 
		e <- expand.grid(1:(segments-1), 1:segments)
		i1 <- apply(e, 1, function(z) z[1] + segments*(z[2] - 1))
		i2 <- i1 + 1
		i3 <- (i1 + segments - 1) %% segments^2 + 1
		i4 <- (i2 + segments - 1) %% segments^2 + 1
		i <- rbind(i1, i2, i4, i3)
		x <- rgl::asEuclidean(t(v))
		ellips <- rgl::qmesh3d(v, i)
		# override settings for 1 df line
		if (df<2) {
			wire <- TRUE
			shade <- FALSE
		}
		if (verbose) print(paste("col:", col, " shade:", shade, " alpha:", alpha, " wire:", wire, sep=" "))
		if(shade) rgl::shade3d(ellips, col=col, alpha=alpha, lit=TRUE)
		if(wire) rgl::wire3d(ellips, col=col, size=lwd, lit=FALSE)
		bbox <- matrix(rgl::par3d("bbox"), nrow=2)
		ranges <- apply(bbox, 2, diff)
		if (!is.null(label) && label !="")
			rgl::texts3d(x[which.max(x[,2]),] + offset*ranges, adj=0, texts=label, color=col, lit=FALSE)
		rownames(bbox) <- c("min", "max")
		return(bbox)
	}
	
	
	#if (!require(car)) stop("car package is required.")
#	if (!require(rgl)) stop("rgl package is required.")    
	if (!requireNamespace("rgl")) stop("rgl package is required.")    
	# avoid deprecated warnings from car
	if (car2 <- packageDescription("car")[["Version"]] >= 2) linear.hypothesis <- linearHypothesis

	type <- match.arg(type)
	size <- match.arg(size)
	fogtype <- match.arg(fogtype)
	bg.col <- match.arg(bg.col)    
	data <- model.frame(mod)
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
	if (verbose) print(manova)    
#	response.names <- rownames(manova$SSPE)
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
	if (length(variables) != 3) {
		extra <- if (length(variables) == 1) 'heplot1d()' else 
				if (length(variables) == 2) 'heplot()' else 'pairs()'
		stop(paste("You may only plot 3 response variables. Use", extra))
	}
	
	if (missing(terms) || (is.logical(terms) && terms)) {
		terms <- manova$terms
		if (!is.null(iterm)) {
			terms <- terms[grep(iterm, terms)]   ## only include those involving iterm
		}
		if (remove.intercept) terms <- terms[terms != "(Intercept)"]
	}
	n.terms <- if (!is.logical(terms)) length(terms) else 0 
	# note: if logical here, necessarily FALSE
	n.hyp <- if (missing(hypotheses)) 0 else length(hypotheses)
	n.ell <- n.terms + n.hyp
	if (n.ell == 0) stop("Nothing to plot.")
	
	E <- SSPE
	p <- nrow(E)
	E <- E[variables, variables]
	Y <- Y[,vars] 
	gmean <- if (missing(data))  c(0,0,0) 
			else colMeans(Y)
	if (missing(xlab)) xlab <- vars[1]
	if (missing(ylab)) ylab <- vars[2]
	if (missing(zlab)) zlab <- vars[3]
	dfe <- manova$error.df
	scale <- 1/dfe 
	E <- E * scale
	radius <- sqrt(3 * qf(level, 3, dfe))
	
	col   <- he.rep(col, n.ell); E.col<- col[length(col)]
	shade <- he.rep(shade, n.ell)
	shade.alpha <- he.rep(shade.alpha, n.ell)
	wire  <- he.rep(wire, n.ell)
	
if (!add){    
		rgl::open3d()
		rgl::view3d(fov=fov)
		rgl::bg3d(color=bg.col, fogtype=fogtype)    
	} 
	
	if (error.ellipsoid) {
		E.ellipsoid <- ellipsoid(gmean, E, radius, col=E.col, label=err.label, 
				shade=shade[[length(shade)]], alpha=shade.alpha[[length(shade.alpha)]],
				wire=wire[[length(wire)]])
		colnames(E.ellipsoid) <- vars
	}       
	term.labels <- if (n.terms == 0) NULL
			else if (!is.logical(term.labels)) term.labels
			else if (term.labels) terms else rep("", n.terms)
	
	H.ellipsoid <- as.list(rep(0, n.ell))
	if (n.terms > 0) for (term in 1:n.terms){
			term.name <- terms[term] 
			H <- manova$SSP[[term.name]]
			H <- H[variables, variables]
			dfh <- manova$df[term.name]
			#          scale <- eval(parse(text=h.scale))
			factor <- if (size == "evidence") lambda.crit(alpha, p, dfh, dfe) else 1  
			H <- H * scale/factor
			if (verbose){
				cat(term.name, " H matrix (", dfh, " df):\n")
				print(H)
			}
			if((!shade[term]) & !wire[term]) 
				warning(paste("shate and wire are both FALSE for ", term), call.=FALSE)
			H.ellipsoid[[term]] <- ellipsoid(gmean, H, radius, col=col[term], label=term.labels[term], 
					df=dfh, shade=shade[term], alpha=shade.alpha[term], wire=wire[term])  
			colnames(H.ellipsoid[[term]]) <- vars
		}
	hyp.labels <- if (n.hyp == 0) NULL
			else if (!is.logical(hyp.labels)) hyp.labels
			else if (hyp.labels) names(hypotheses) else rep("", n.hyp)  
	if (n.hyp > 0) for (hyp in 1:n.hyp){
			lh <- linearHypothesis(mod, hypotheses[[hyp]])
			H <- lh$SSPH[variables, variables]
			dfh <- lh$df
			factor <- if (size == "evidence") lambda.crit(alpha, p, dfh, dfe) else 1  
			H <- H * scale/factor
			if (verbose){
				cat("\n\n Linear hypothesis: ", names(hypotheses)[[hyp]], "\n") 
				print(lh)
			}
			term <- n.terms + hyp
			H.ellipsoid[[term]] <- ellipsoid(gmean, H, radius, col=col[term], label=hyp.labels[hyp],
					df=dfh, shade=shade[term], alpha=shade.alpha[term], wire=wire[term])
		}         
	ranges <- apply(matrix(rgl::par3d("bbox"), nrow=2), 2, diff)
#   if (grand.mean) ellipsoid(gmean, diag((ranges/40)^2), col="black", wire=FALSE, alpha=0.8) # centre dot    
	# better: use a centered 3D cross here
	if (grand.mean) cross3d(gmean, (ranges/25), col="black", lwd=2) # centre cross            
	
#browser()	
	## BUG fixed here:  should only label the means for factors included in terms
	if ((!is.logical(factor.means)) || factor.means){
		factors <- data[, sapply(data, is.factor), drop=FALSE]
		factor.names <- colnames(factors) 
		if (is.null(iterm)) factor.names <- factor.names[factor.names %in% terms]
		if (!is.logical(factor.means)){
			which <- match(factor.means, factor.names)
			check <- is.na(which)
			if (any(check)) stop(paste(factor.means[check], collapse=", "), 
						" not among factors.")
			factors <- factors[, which, drop=FALSE]
		}
		else factors <- factors[, factor.names, drop=FALSE]    
#        for (fac in factors){
#            means <- aggregate(Y, list(fac), mean)
		if (ncol(factors)) for (j in 1:ncol(factors)){
				means <- aggregate(Y, list(factors[,j]), mean)
				# color the points the same as the ellipse for the term
				loc <- match(factor.names[j], terms, nomatch=0)
				pcol <- if (loc>0) col[loc] else "black"
				for (m in 1:nrow(means)) {
					ellipsoid(unlist(means[m, 2:4]), diag((ranges/100))^2, col=pcol, wire=FALSE, alpha=0.8)
#            		points3d(unlist(means[m, 2:4]), size=3, color=pcol)
#					spheres3d(unlist(means[m, 2:4]), radius=diag((ranges/30))^2, color=pcol)
				}
				rgl::texts3d(means[,2:4] + matrix(offset*ranges, nrow(means), 3, byrow=TRUE), 
						texts=as.character(means[,1]), color=pcol, adj=0)
			}
	}
	
	# handle xlim, ylim, zlim
	## enforce that the specified limits are at least as large as the bbox
	if (!missing(xlim) | !missing(ylim) | !missing(zlim)) {
		bbox <- matrix(rgl::par3d("bbox"),3,2,byrow=TRUE)
		xlim <- if(missing(xlim)) bbox[1,] else c(min(xlim[1],bbox[1,1]), max(xlim[2],bbox[1,2]))
		ylim <- if(missing(ylim)) bbox[2,] else c(min(ylim[1],bbox[2,1]), max(ylim[2],bbox[2,2]))
		zlim <- if(missing(zlim)) bbox[3,] else c(min(zlim[1],bbox[3,1]), max(zlim[2],bbox[3,2]))
		rgl::decorate3d(xlim=xlim, ylim=ylim, zlim=zlim, box=FALSE, axes=FALSE, xlab=NULL, ylab=NULL, zlab=NULL, top=FALSE)
	}
	
	if (add) rgl::rgl.pop(id=savedvars$.frame)
	frame <- rgl::axis3d("x-", color="black")
	frame <- c(frame, rgl::mtext3d(xlab, "x-", color="black", line=1.5))
	frame <- c(frame, rgl::axis3d("y-", col="black"))
	frame <- c(frame, rgl::mtext3d(ylab, "y-", color="black", line=1.5))
	frame <- c(frame, rgl::axis3d("z-", col="black"))
	frame <- c(frame, rgl::mtext3d(zlab, "z-", color="black", line=1.5))
	frame <- c(frame, rgl::box3d(col="black"))
	assign(".frame", frame, envir=savedvars)
	#   savedvars$.frame <- frame

	rgl::aspect3d(x=1, y=1, z=1)
	
	names(H.ellipsoid) <- c(if (n.terms > 0) term.labels, if (n.hyp > 0) hyp.labels)
	result <- if(error.ellipsoid) list(H=H.ellipsoid, E=E.ellipsoid, center=gmean, radius=radius) 
			else list(H=H.ellipsoid, center=gmean, radius=radius)
	class(result) <- "heplot3d"
	invisible(result)
	
}

.frame <- NULL   # avoid warning

