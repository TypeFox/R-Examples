#' Sequential loglinear models for an n-way table

#' This function takes an n-way contingency table and fits a series of sequential
#' models to the 1-, 2-, ... n-way marginal tables, corresponding to a variety of
#' types of loglinear models.


#' @param x 	a contingency table in array form, with optional category labels specified in the dimnames(x) attribute,
#'            or else a data.frame in frequency form, with the frequency variable names "Freq".
#' @param type type of sequential model to fit
#' @param marginals which marginals to fit?
#' @param vorder  order of variables
#' @param k    indices of conditioning variable(s) for "joint", "conditional" or order for "markov"
#' @param prefix
#' @param fitted keep fitted values?


seq_loglm <- function(
	x,
	type = c("joint", "conditional", "mutual", "markov", "saturated"),
	marginals = 1:nf,  # which marginals to fit?
	vorder = 1:nf,     # order of variables in the sequential models
	k = NULL,          # conditioning variable(s) for "joint", "conditional" or order for "markov"
	prefix = 'model',
	fitted = TRUE,     # keep fitted values?
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
#  models <- as.list(rep(NULL, length(marginals))) 
  models <- list()
  for (i in marginals) {
		mtab <- margin.table(x, 1:i)
		if (i==1) {
			# KLUDGE: use loglin, but try to make it look like a loglm object
			mod <- loglin(mtab, margin=NULL, print=FALSE)
		  mod$model.string = paste("=", factors[1])
		  mod$margin <- list(factors[1])
#		  mod$margin <- names(dimnames(mtab))
#		  names(mod$margin) <- factors[1]
      if (fitted) {
        fit <- mtab
  		  fit[] <- (sum(mtab) / length(mtab))
  		  mod$fitted <- fit
		  }
		  mod$nobs <- length(mtab)
		  mod$frequencies <- mtab
		  mod$deviance <- mod$lrt
		  class(mod) <- c("loglin", "loglm")
		  }
		else {
  		expected <- switch(type,
  			'conditional' = conditional(i, mtab, with=if(is.null(k)) i else k),
  			'joint' = joint(i, mtab, with=if(is.null(k)) i else k),
  			'mutual' = mutual(i, mtab),
  			'markov' = markov(i, mtab, order=if(is.null(k)) 1 else k),
  			'saturated' = saturated(i, mtab)
  			)

  		form <- loglin2formula(expected)
#  		mod <- loglm(formula=form, data=mtab, fitted=TRUE)
      mod <- eval(bquote(MASS::loglm(.(form), data=mtab, fitted=fitted)))
  		mod$model.string <- loglin2string(expected, brackets=if (i<nf) '()' else '[]')
  		
		}
#  	cat(i, "  model.string: ", mod$model.string, "\n")
#  	cat("model:\n"); print(mod)
  	models[[i]] <- mod
	}
	names(models) <- paste(prefix, '.', indices, sep='')
	models <- models[marginals]
	class(models) <- "loglmlist"
	invisible(models)
}

