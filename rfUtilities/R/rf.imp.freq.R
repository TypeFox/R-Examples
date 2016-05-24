#' @title Random Forest variable selection frequency
#' @description Evaluates the frequency that an independent variables are selected greater-than/equal-to defined importance threshold
#'    
#' @param x random forest object
#' @param p Threshold of row standardized importance values
#' @param plot Plot frequencies (TRUE/FALSE)
#' 
#' @return A list class object with the following components:
#'  frequency:
#'    vars - names of independent variables used in model
#'    global - if a variable greater-than/equal-to importance threshold, else NA 
#'    column for each class where greater-than/equal-to importance threshold, else NA    
#'    var.freq - frequency a variable is selected for global and local importance >= importance threshold
#'
#'  importance: 
#'    Standardized importance matrix from randomForest model 
#'
#' @note
#'   Evaluates the number of times a variable is selected greater-than/equal-to defined threshold (p) for the global and local (class level) importances. This allow one to evaluate if a given variable is important to the overall model or specific classes.  
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @examples  
#'  require(randomForest)
#'  data(iris)
#'  iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE)
#'  rf.imp.freq(iris.rf, p = 0.30)
#'
#' @export
rf.imp.freq <- function(x, p = 0.60, plot = TRUE) {
  if (!inherits(x, "randomForest")) stop("x is not randomForest class object")
  if (x$type == "classification" | x$type == "unsupervised") {
   if (is.null(x$importanceSD) == TRUE | "MeanDecreaseAccuracy" %in% 
      names(as.data.frame(x$importance)) == FALSE)
      stop("randomForest object does not contain importance, please run with importance=TRUE")  
	} else {
	  stop("Does not support regression")
    }	
  imp.measure <- "MeanDecreaseAccuracy"  
  importance <- x$importance
    for(i in 1:ncol(importance)) { 
	  importance[,i] <- importance[,i] / max(importance[,i])
    }
  vars <- row.names(importance)
  imp <- row.names(importance[importance[,imp.measure] > p ,]) 
  tmp.index <- which(is.na(match(vars,imp)))    
  tmp.vars <- vars
  tmp.vars[tmp.index] <- NA		 
  vars.df <- data.frame(vars=vars, global=tmp.vars)
    for(j in 1:length(x$classes)) { 
       x.class <- x$classes[j]  
       tmp <- row.names(importance[importance[,x.class] > p ,]) 
       tmp.index <- which(is.na(match(vars,tmp))) 
       tmp.vars <- vars
       tmp.vars[tmp.index] <- NA		 
       vars.df <- data.frame(vars.df, tmp.vars)
         names(vars.df)[j+2] <- x.class
      }
	vars.df <- data.frame(vars.df, var.freq=apply(vars.df[,2:ncol(vars.df)], MARGIN=1, 
                          FUN=function(x) { length(x[!is.na(x)]) }  ) )  
	if(plot == TRUE) {
	  graphics::barplot(vars.df$var.freq, names.arg=vars.df$vars, 
	              main="Frequency of variable importance", las=2)
	}
  return( list( frequency = vars.df, importance = importance) )
  }
