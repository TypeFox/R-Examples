

# Create longitudinal dataset.
# Instantiate new ggobi with a longitudinal data set.
#
# This function allows you to load longitudinal data in
# to GGobi and display it as a line plot.  This is achieved
# by creating edges between adjacent time points, for a given
# id variable.
#
# For best viewing, we recommend that you turn the show points off
# in the options menu.  When brushing, you may also want to use
# categorical brushing on the id variable, so that the entire
# series is selected for an observation.
#
# @arguments data frame
# @arguments time variable
# @arguments id variable
# @arguments ggobi instance, if you don't want to create a new one
# @keyword dynamic
#X data(Oxboys, package="nlme")
#X ggobi_longitudinal(Oxboys, Occasion, Subject)
#X ggobi_longitudinal(stormtracks, seasday, id)
#X ggobi_longitudinal(data.frame(x=1:100, y=sin(1:100)))
ggobi_longitudinal <- function(data, time=1:rows, id=rep(1, rows), g = NULL) {
	name <- deparse(substitute(data))
	rows <- nrow(data)
	time <- eval(substitute(time), data)
	obsUnit <- eval(substitute(id), data)

	or <- order(obsUnit, time)
	tmp <- data[or, ]
	if (is.null(g)) {
	  g <- ggobi(tmp, name=name)
	} else {
	  g[name] <- tmp
	}

	edges <- cbind(rownames(tmp[-nrow(tmp), ]), rownames(tmp[-1, ]))
	matching <- obsUnit[or][-1] == obsUnit[or][-nrow(tmp)]
	edges[!matching, ] <- NA

	edges <- edges[complete.cases(edges),]

	d <- data.frame(tmp[edges[,1], sapply(tmp, is.factor), drop=FALSE])
	rownames(d) <- paste(name, 1:nrow(d), sep="")
	g[paste(name, "edges", sep="-")] <- d
	edges(g[paste(name, "edges", sep="-")]) <- edges

	if (!is.null(ggobi)) {
    d <- displays(g)[[1]]
    edges(d) <- g[2]
  }

	invisible(g)
}

# Create parallel coordinates plot.
# Mock up a pcp plot using points and edges.
#
# Experimental and may suggest ways to reduce PCP code in GGobi
#
# @keyword internal
ggobi_pcp <- function(data, type="range") {
  if (!require(reshape)) stop("Must have reshape package installed")
  data$CASEID <- factor(1:nrow(data))
  datam <- melt(rescaler(data, type=type), id = "CASEID")

  g <- ggobi_longitudinal(datam, id=CASEID)
  g$raw <- data

  d <- displays(g)[[1]]
  edges(d) <- g[1]
  variables(d) <- list(X="variable", Y="value")

  invisible(g)
}
globalVariables("CASEID")
