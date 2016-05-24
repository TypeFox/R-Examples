#' enter
#'
#' Helper Function for xpssRegression
#'
#' Regression Method. All stated variables are used at once in the Regression Modell.
#' This function only works inside xpssRegression \code{\link{xpssRegression}} 
#' @usage enter(...) 
#' @param ... one or more Variables specified as characters
#' @return returns a list of class summary.lm
#' @author Martin Schneider
#' @seealso xpssRegression \code{\link{xpssRegression}}
#' @examples \dontrun{
#' enter("var1","var2","var3")
#' }
#' @keywords internal
#' @export
enter <- function(...){
	options(warnings = -1)
	variables <- c(...)
	
	if(is.null(variables)){
		variables <- get("variables",parent.frame())
	}
	
	if(length(variables) == 1 & variables == "ALL"){
		
	}
	
	dependent <- get("dependent", parent.frame())
	
	x <- get("x", parent.frame())
	origin <- get("origin", parent.frame())
	
	#formel bilden
	formel <- paste0(dependent, "~", paste0(variables, collapse = "+"))
	
	if(origin)formel <- paste0(formel,"+0")
	
	#eigentliche Analyse
	reg <- lm(as.formula(formel), data = x)		
	reg$call <- paste0("lm(",formel,", data = x")
	options(warnings = 0)
	return(reg)
}


#' Calculates a linear Regression
#'
#' R Implementation of the SPSS \code{REGRESSION} Function.xpssRegression calculates linear regressions with associated statistics and plots. 
#'
#' Implementation from SPSS Regression in R
#'
#' @usage xpssRegression(x, variables = NULL, dependent, method = enter(), 
#' missing = "listwise", statistics = c("COEFF", "OUTS", "R", "ANOVA"), origin = FALSE)
#' @param x a (non-empty) data.frame, data.table object or input data of class "xpssFrame". 
#' @param variables vector with independent Variables as characters
#' @param dependent dependent Variable as character
#' @param method regression method to apply. Currently only \code{enter()} is implemented
#' @param missing Method to handle missing values. Currently only listwise is implemented
#' @param statistics character vector, which statistics should be in the output. Currently the function will return standard \code{\link{lm}} output. Only output for ANOVA is optional
#' @param origin Should the constant be compressed? 
#' @return returns a list of lists, each applied method is its own list with regression, anova and beta-coeffizients as elements
#' @author Martin Schneider
#' @examples
#' 
#' data(fromXPSS)
#' 
#' xpssRegression(x = fromXPSS,
#'  		variables = c("V7_1","V7_2"),
#'  		dependent = "V5",
#'  		method = list(enter()))
#' 
#' @export

xpssRegression <- function(x,
		variables = NULL,
		dependent,
		method = enter(),
		missing = "listwise",
		statistics = c("COEFF", "OUTS", "R", "ANOVA"),
		origin = FALSE
)
{
  
	
  ####################################################################
  ####################################################################
  
  functiontype <- "AN"
  x <- applyMetaCheck(x)
  
  ####################################################################
  ####################################################################
  ####################################################################
  
  
	options(warn = -1)
	
	method <- substitute(method)
	method <- eval(list(method)[[1]])
	
	force(x)
  
  #include
  if(missing == "include") {
    x <- computeValue(x, variables)
  }
	
	#listwise
	if(missing == "listwise") {
		x <- na.omit(x)
	}
  if(missing == "pairwise") {
    
  }
	
	result <- list()
	for(i in 1:length(method)){
		method_i <- method[[i]]
		result[[i]] <- list()
		
		names(result)[i] <- paste(names(method_i$model)[1],"~", paste(names(method_i$model)[-1],collapse = "+"))
		
		if("ANOVA"%in%statistics)result[[i]]$anova <- anova(eval(method_i))
				
		result[[i]]$regression <- summary(eval(method_i))
		
		#beta
		b <- summary(method_i)$coef[-1, 1]
		sx <- sapply(method_i$model[-1], sd)
		sy <- sapply(method_i$model[1], sd)
		result[[i]]$beta <- b * sx/sy
	}	
	options(warn = 0)
	return(result)
}

print.xpssRegression <- function(result){
	cat(summary)
}


