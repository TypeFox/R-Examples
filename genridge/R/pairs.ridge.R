pairs.ridge <-
function(x, variables, radius=1, lwd=1, lty=1,
		col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
		center.pch = 16, center.cex=1.25, digits=getOption("digits") - 3,
		diag.cex= 2, diag.panel = panel.label,
		fill=FALSE, fill.alpha=0.3, ...) {

	panel.label <- function(x, ...) {
		op <- par(xpd=TRUE)
		on.exit(par(op))
		plot(c(min, max),c(min, max), type="n", axes=FALSE)
		text(0.5, 0.5, vars[i], cex=diag.cex)
		text(1, 0, signif(range[1, i], digits=digits), adj=c(1, 0))
		text(0, 1, signif(range[2, i], digits=digits), adj=c(0, 1)) 
		box()
	}	

#	why doesn't this work??
	panel.barplot <- function(x, ...) {
		barplot(x$coef[,i], axes=FALSE, col=col, ...)
		box()
	}

	vars <- dimnames(x$coef)[[2]]
  if (!missing(variables)){
      if (is.numeric(variables)) {
          vars <- vars[variables]
          if (any(is.na(vars))) stop("Bad variable selection.")
          }
      else {
          check <- !(variables %in% vars)
          if (any(check)) stop(paste("The following", 
              if (sum(check) > 1) "variables are" else "variable is",
              "not in the model:", paste(variables[check], collapse=", ")))
          vars <- variables
          }
      }
  else variables <- vars
	nvar <- length(vars)

  range <- apply(x$coef, 2, range)
	min=0
	max=1	
  old.par <- par(mfrow=c(nvar, nvar), mar=rep(0,4))
  on.exit(par(old.par))

#	browser()
	for (i in 1:nvar){
    for (j in 1:nvar){
      if (i == j)
				diag.panel(x)
      else {
        plot.ridge(x, variables=c(vars[j], vars[i]), radius=radius,
			labels=NULL,
        	col=col, lwd=lwd, lty=lty, center.cex=center.cex,
        	axes=FALSE,
          fill=fill, fill.alpha=fill.alpha, ...)
        box()
            }
        }
    }

}

