# last modified 21 December 2006 by J. Fox
# last modified 15 April 2009 by M. Friendly 
#  -- Fixed numerous warnings resulting from axes=FALSE
#  -- prepare to generalize diagonal panel
# last modified 2 Feb 2010 by M. Friendly
#  -- added code for repeated measures designs
# last modified 13 Oct 2011 by M. Friendly
#  -- added var.labels 


`pairs.mlm` <-
function(x, variables, var.labels, var.cex = 2,
    type=c("II", "III", "2", "3"),
	idata=NULL,
	idesign=NULL,
	icontrasts=NULL,
	imatrix=NULL,
	iterm=NULL,
	manova,        # an optional Anova.mlm object
	offset.axes=0.05, 
	digits=getOption("digits") - 1,
	fill=FALSE,         ## whether to draw filled ellipses (vectorized)
	fill.alpha=0.3,     ## alpha transparency for filled ellipses
	...){

#	manova <- Anova(x, type)
	if (missing(manova)) {
		type <- match.arg(type)
		if (is.null(imatrix)) {
			manova <- Anova(x, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts)
		}
		else {
			if (packageDescription("car")[["Version"]] >= 2)
				manova <- Anova(x, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)
			else stop("imatrix argument requires car 2.0-0 or later")
		} 
	}   
	
	data <- model.frame(x)
#	Y <- model.response(model.frame(x))
	if (is.null(idata) && is.null(imatrix)) {
		Y <- model.response(data) 
#		SSPE <- manova$SSPE
	} 
	else {
		if (is.null(iterm)) stop("Must specify a within-S iterm for repeated measures designs" )
		### FIXME::car -- workaround for car::Anova.mlm bug: no names assigned to $P component
		if (is.null(names(manova$P))) names(manova$P) <- names(manova$SSPE)
		Y <- model.response(data) %*% manova$P[[iterm]]
#		SSPE <- manova$SSPE[[iterm]]
	}   
	
	vars <- colnames(Y)
  if (!missing(variables)){
      if (is.numeric(variables)) {
          vars <- vars[variables]
          if (any(is.na(vars))) stop("Bad response variable selection.")
          }
      else {
          check <- !(variables %in% vars)
          if (any(check)) stop(paste("The following", 
              if (sum(check) > 1) "variables are" else "variable is",
              "not in the model:", paste(variables[check], collapse=", ")))
          vars <- variables
          }
      }

	if(missing(var.labels)) var.labels <- vars
	else {
		if (length(var.labels) < length(vars)) stop("Too few var.labels supplied")
	}
	
    n.resp <- length(vars)
    if (n.resp < 3) stop("Fewer than 3 response variables.")
    range <- apply(Y, 2, range)
    min <- - offset.axes
    max <- 1 + offset.axes
    old.par <- par(mfrow=c(n.resp, n.resp), mar=rep(0,4))
    on.exit(par(old.par))

	panel.label <- function(x, ...) {
		plot(c(min, max),c(min, max), type="n", axes=FALSE)
		text(0.5, 0.5, var.labels[i], cex=var.cex)
		text(1, 0, signif(range[1, i], digits=digits), adj=c(1, 0))
		text(0, 1, signif(range[2, i], digits=digits), adj=c(0, 1)) 
		box()
	}	
	for (i in 1:n.resp){
	  for (j in 1:n.resp){
	    if (i == j){
	      panel.label()
	    }
	    else {
	      heplot(x, variables=c(vars[j], vars[i]), manova=manova, axes=FALSE,
	             idata=idata, idesign=idesign, imatrix=imatrix, iterm=iterm,
	             offset.axes=offset.axes, fill=fill, fill.alpha=fill.alpha, ...)
	      box()
	    }
	  }
	}
}

