#' Sequential Mosaics and Strucplots for an N-way Table


#' This function takes an n-way contingency table and plots mosaics for series of sequential
#' models to the 1-, 2-, ... n-way marginal tables, corresponding to a variety of
#' types of loglinear models.

#' @param x 	a contingency table in array form, with optional category labels specified in the dimnames(x) attribute,
#'            or else a data.frame in frequency form, with the frequency variable names "Freq".
#' @param panel  panel function
#' @param type type of sequential model to fit
#' @param plots which marginals to plot?
#' @param vorder  order of variables
#' @param k    indices of conditioning variable(s) for "joint", "conditional" or order for "markov"
#' @export

seq_mosaic <- function(
	x,
	panel = mosaic, 
	type = c("joint", "conditional", "mutual", "markov", "saturated"),
	plots = 1:nf,      # which plots to produce?
	vorder = 1:nf,      # order of variables in the sequential plots
	k = NULL,          # conditioning variable(s) for "joint", "conditional" or order for "markov"
	...
	)
{

  if (inherits(x, "data.frame") && "Freq" %in% colnames(x)) {
    x <- xtabs(Freq ~ ., data=x)
  }
  if (!inherits(x, c("table", "array"))) stop("not an xtabs, table, array or data.frame with a 'Freq' variable")
  
	nf <- length(dim(x))
	x <- aperm(x, vorder)
	factors <- names(dimnames(x))
	indices <- 1:nf
	
  type = match.arg(type)
	for (i in plots) {
		mtab <- margin.table(x, 1:i)
		df <- NULL
		if (i==1) {
		  expected <- mtab
		  expected[] <- sum(mtab) / length(mtab)
		  df <- length(mtab)-1
		  model.string = paste("=", factors[1])
		  }
		else {
  		expected <- switch(type,
  			'conditional' = conditional(i, mtab, with=if(is.null(k)) i else k),
  			'joint' = joint(i, mtab, with=if(is.null(k)) i else k),
  			'mutual' = mutual(i, mtab),
  			'markov' = markov(i, mtab, order=if(is.null(k)) 1 else k),
  			'saturated' = saturated(i, mtab)
  			)
    model.string <- loglin2string(expected, brackets=if (i<nf) '()' else '[]')
  	}
		panel(mtab, expected=expected, df=df, main=model.string, ...)
	}
}

