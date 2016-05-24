#' @export
LSresids <- function(object, ...){
  UseMethod("LSresids", object)
}

#' @export
#' @rdname LSresids.mer
#' @method LSresids default
LSresids.default <- function(object, ...){
  stop(paste("there is no LSresids() method for objects of class",
             paste(class(object), collapse=", ")))
}

#' Calculating least squares residuals
#'
#' This function calculates least squares (LS) residuals
#' found by fitting separate LS regression models to each case.
#' For examples see the documentation for \code{HLMresid}.
#'
#' @export
#' @method LSresids mer
#' @S3method LSresids mer
#' @aliases LSresids
#' @param object an object of class \code{mer} or \code{lmerMod}.
#' @param level which residuals should be extracted: 1 for case-level
#' residuals or the name of a grouping factor (as defined in \code{flist} of the 
#' \code{mer} object) for between-group residuals.
#' @param sim optional argument giving the data frame used for LS residuals. This
#'  is used mainly when dealing with simulations.
#' @param standardize if \code{TRUE} the standardized level-1
#' residuals will also be returned (if \code{level = 1}); if \code{"semi"} then
#' the semi-standardized level-1 residuals will be returned.
#' @param ... do not use
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @references 
#' Hilden-Minton, J. (1995) Multilevel diagnostics for mixed and hierarchical 
#' linear models. University of California Los Angeles.
#' @export
#' @seealso \code{\link{HLMresid}}
#' @keywords models regression
LSresids.mer <- function(object, level, sim = NULL, standardize = FALSE, ...){
  if(!object@dims["LMM"]){
    stop("LSresids is currently not implemented for GLMMs or NLMMs.")
  }
  if(!level %in% c(1, names(object@flist))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  if(object@dims[["nest"]] == 0) {
    stop("LSresids has not yet been implemented for models with 
          crossed random effects")
  }
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi' .")
  }
  
  LS.resid <- NULL # Make codetools happy
  
  fixed <- as.character( fixform( formula(object) ) )
	
	data <- object@frame
	if(!is.null(sim)){data[,fixed[2]] <- sim}
	
	if(level == 1){
		# fitting a separate LS regression model to each group
		form <- paste(fixed[2], fixed[1], fixed[3], "|", names(object@flist)[1])
	
		ls.models <- adjust_lmList(object = formula(form), data = data)
		
		ls.residuals <- lapply(ls.models, resid)
		ls.fitted <- lapply(ls.models, fitted)
	
		# creating a data frame of the residuals, fitted values, and model frames
		ls.data <- lapply(ls.models, model.frame)
		res.data <- do.call('rbind', ls.data)
		
		row.order <- unlist(lapply(ls.data, function(x) row.names(x)))

# 		row.order <- as.numeric(unlist(lapply(ls.data, function(x) row.names(x))))
		
		return.df <- data.frame(LS.resid = unlist(ls.residuals), 
                            fitted = unlist(ls.fitted))
	
		if(!is.null(standardize) && standardize == "semi"){
		  ls.influence <- lapply(ls.models, lm.influence)
		  ls.hat <- lapply(ls.influence, function(x) x$hat)
      
		  h <- unlist(ls.hat)
			semi.std.resid  <- with(return.df, LS.resid / sqrt(1 - h))
			semi.std.resid[is.infinite(semi.std.resid)] <- NaN
	
			return.df <- cbind(return.df, semi.std.resid = semi.std.resid)
		}
    
		if(!is.null(standardize) && standardize == TRUE){
		  ls.rstandard <- unlist(lapply(ls.models, rstandard))
		  ls.rstandard[is.infinite(ls.rstandard)] <- NaN
      
      return.df <- cbind(return.df, std.resid = ls.rstandard)
		}
		
# 		return.df <- return.df[order(row.order),]
		return.df <- cbind(data, return.df)
    rownames(return.df) <- row.order

		return(return.df)
	}
	
	if(level != 1){
    n.ranefs <- length(names(object@flist))
		ranef_names <- names( lme4::ranef(object)[[level]] )
		
		form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
    
    if( level == names(object@flist)[n.ranefs] ) { # For highest level unit
      gform <- fixform( formula(object) )
      global.model <- lm(formula = formula(gform), data = data)
      
      ls.models <- adjust_lmList(object = formula(form), data = data)
      ls.resid <- coef(ls.models)[,ranef_names] - coef(global.model)[ranef_names]
      
    } else{ # For 'intermediate' level unit
      higher.level <- names(object@flist)[which(names(object@flist) == level) + 1]
      gform <- paste(fixed[2], fixed[1], fixed[3], "|", higher.level)
      global.model <- adjust_lmList(object = formula(gform), data = data)
      
      ls.models <- adjust_lmList(object = formula(form), data = data)
      
      # matching lower level units to higer level units
      units.mat <- unique(object@flist)
      n.higher.level <- table(units.mat[higher.level])
      
      global.coefs <- rep(unlist(coef(global.model)[ranef_names]), 
                          times = n.higher.level)
      ls.resid <- coef(ls.models)[,ranef_names] - global.coefs
      
    }
		
		ls.resid <- as.data.frame(ls.resid)
		colnames(ls.resid) <- ranef_names

		return(ls.resid)
	}
}

# 'fixform' is a copy of the 'nobars' function in the lme4 package, 
# renamed so it doesn't cause any conflicts.  This should not be exported.
fixform <- function (term) 
{
  if (!("|" %in% all.names(term))) 
    return(term)
  if (is.call(term) && term[[1]] == as.name("|")) 
    return(NULL)
  if (length(term) == 2) {
    nb <- fixform(term[[2]])
    if (is.null(nb)) 
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- fixform(term[[2]])
  nb3 <- fixform(term[[3]])
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}



#' @export
#' @rdname LSresids.mer
#' @method LSresids lmerMod
#' @S3method LSresids lmerMod
LSresids.lmerMod <- function(object, level, sim = NULL, standardize = FALSE, ...){
  if(!lme4::isLMM(object)){
    stop("LSresids is currently not implemented for GLMMs or NLMMs.")
  }
  if(!level %in% c(1, names(object@flist))) {
    stop("level can only be 1 or a grouping factor from the fitted model.")
  }
  if(!isNestedModel(object)) {
    stop("LSresids has not yet been implemented for models with 
         crossed random effects")
  }
  if(!is.null(standardize) && !standardize %in% c(FALSE, TRUE, "semi")) {
    stop("standardize can only be specified to be logical or 'semi' .")
  }
  
  LS.resid <- NULL # Make codetools happy
  
  fixed <- as.character( fixform( formula(object) ) )
  
  data <- object@frame
  if(!is.null(sim)){data[,fixed[2]] <- sim}
  
  if(level == 1){
    # fitting a separate LS regression model to each group
    form <- paste(fixed[2], fixed[1], fixed[3], "|", names(object@flist)[1])
    
    ls.models <- adjust_lmList(object = formula(form), data = data)
    
    ls.residuals <- lapply(ls.models, resid)
    ls.fitted <- lapply(ls.models, fitted)
    
    # creating a data frame of the residuals, fitted values, and model frames
    ls.data <- lapply(ls.models, model.frame)
    res.data <- do.call('rbind', ls.data)
    
    row.order <- unlist(lapply(ls.data, function(x) row.names(x)))
    
    return.df <- data.frame(LS.resid = unlist(ls.residuals), 
                            fitted = unlist(ls.fitted))
    
    if(!is.null(standardize) && standardize == "semi"){
      ls.influence <- lapply(ls.models, lm.influence)
      ls.hat <- lapply(ls.influence, function(x) x$hat)
      
      h <- unlist(ls.hat)
      semi.std.resid  <- with(return.df, LS.resid / sqrt(1 - h))
      semi.std.resid[is.infinite(semi.std.resid)] <- NaN
      
      return.df <- cbind(return.df, semi.std.resid = semi.std.resid)
    }
    
    if(!is.null(standardize) && standardize == TRUE){
      ls.rstandard <- unlist(lapply(ls.models, rstandard))
      ls.rstandard[is.infinite(ls.rstandard)] <- NaN
      
      return.df <- cbind(return.df, std.resid = ls.rstandard)
    }
    
    return.df <- cbind(res.data, return.df)
    rownames(return.df) <- row.order
    
    return(return.df)
  }
  
  if(level != 1){
    n.ranefs <- length(names(object@flist))
    ranef_names <- names( lme4::ranef(object)[[level]] )
    
    form <- paste(fixed[2], fixed[1], fixed[3], "|", level)
    
    if( level == names(object@flist)[n.ranefs] ) { # For highest level unit
      gform <- fixform( formula(object) )
      global.model <- lm(formula = formula(gform), data = data)
      
      ls.models <- adjust_lmList(object = formula(form), data = data)
      ls.resid <- coef(ls.models)[,ranef_names] - coef(global.model)[ranef_names]
      
    } else{ # For 'intermediate' level unit
      higher.level <- names(object@flist)[which(names(object@flist) == level) + 1]
      gform <- paste(fixed[2], fixed[1], fixed[3], "|", higher.level)
      global.model <- adjust_lmList(object = formula(gform), data = data)
      
      ls.models <- adjust_lmList(object = formula(form), data = data)
      
      # matching lower level units to higer level units
      units.mat <- unique(object@flist)
      n.higher.level <- table(units.mat[higher.level])
      
      global.coefs <- rep(unlist(coef(global.model)[ranef_names]), 
                          times = n.higher.level)
      ls.resid <- coef(ls.models)[,ranef_names] - global.coefs
      
    }
    
    ls.resid <- as.data.frame(ls.resid)
    colnames(ls.resid) <- ranef_names
    
    return(ls.resid)
  }
  }
