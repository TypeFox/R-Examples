### Various contributed functions


##################################################################
################         sort.data.frame         #################
##################################################################


#'Sort a data.frame
#'
#'Sorts a data frame by one or more variables
#'
#'
#'@param x Data.frame to sort
#'@param formula Formula by which to sort the data.frame (e.g. ~group1+group2
#'sorts first by group1 then by group2)
#'@param decreasing Ignored.  Exists for compatibility with generic S3 method.
#'@param \dots Used to pass ,drop=FALSE to [
#'@return Returns a sorted data.frame
#'@note Modifications by Ari Friedman and Roman Lustrik Original Author: Kevin
#'Wright http://tolstoy.newcastle.edu.au/R/help/04/09/4300.html Some ideas from
#'Andy Liaw http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html Use + for
#'ascending, - for decending. Sorting is left to right in the formula
#'
#'If you are Kevin Wright, please contact me.  I have attempted to reach you by
#'every means thinkable, to no avail.  My assumption is that this is in the
#'public domain since you posted it for others to use, but please tell me if
#'that is not the case.
#'@author Kevin Wright, with generic compatibility by Ari B. Friedman
#'@seealso \link[plyr]{arrange}
#'@export sort.data.frame
#'@method sort data.frame
#'@S3method sort data.frame
#'@examples
#'
#'library(datasets)
#'sort.data.frame(ChickWeight,formula=~weight+Time)
#'
#'mydf <- data.frame(col1 = runif(10))
#'rownames(mydf) <- paste("x", 1:10, sep = "")
#'sort(mydf, f = ~col1) # drops a dimension
#'sort(mydf, f = ~col1, drop = FALSE) # does not drop a dimension (returns a data.frame)
#'
#'
sort.data.frame <- function(x, decreasing = NULL, formula, ...) {
	# Author: Kevin Wright
	# http://tolstoy.newcastle.edu.au/R/help/04/09/4300.html
	# Some ideas from Andy Liaw
	# http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html
	# Use + for ascending, - for decending.
	# Sorting is left to right in the formula
	# Useage is either of the following:
	# sort.data.frame(~Block-Variety,Oats)
	# sort.data.frame(Oats,~-Variety+Block)
	
	# If dat is the formula, then switch form and dat
	if(inherits(x,"formula")){
		f=x
		dat=formula
		formula=f
	}
	if(sides(formula)!=1) {
		stop("Formula must be one-sided.")
	}
	# Make the formula into character and remove spaces
	formc <- as.character(formula[2])
	formc <- gsub(" ","",formc)
	# If the first character is not + or -, add +
	if(!is.element(substring(formc,1,1),c("+","-"))) {
		formc <- paste("+",formc,sep="")
	}
	# Extract the variables from the formula
	vars <- unlist(strsplit(formc, "[\\+\\-]"))
	vars <- vars[vars!=""] # Remove spurious "" terms
	# Build a list of arguments to pass to "order" function
	calllist <- list()
	pos=1 # Position of + or -
	for( i in seq(length(vars)) ) {
		varsign <- substring(formc,pos,pos)
		pos <- pos + 1 + nchar(vars[i])
		if( is.factor( x[,vars[i]] ) ) {
			if(varsign=="-")
				calllist[[i]] <- -rank(x[,vars[i]])
			else
				calllist[[i]] <- rank(x[,vars[i]])
		}
		else {
			if( varsign == "-" )
				calllist[[i]] <- -x[,vars[i]]
			else
				calllist[[i]] <- x[,vars[i]]
		}
	  }
	  x[do.call("order",calllist),,...]
} 




##################################################################
################             xtablelm            #################
##################################################################

# Contributed by Mark W. Clements

##### This function produces the output of an lm object as it appears
##### in the R console when you type summary(lmobject)

##### Inputs:
##### lm.obect - is the name of your linear model object that you want to make a summary table for.
##### titref - the label name of the equation you made in Latex to cross reference
##### labname - the label name you want for this table
##### extracaption - adds whatever text string you pass to the title of the table.



#'Produces the output of an lm object as it appears in the R console when you
#'type summary(lmobject)
#'
#'Produces the output of an lm object as it appears in the R console when you
#'type summary(lmobject)
#'
#'
#'@param lm.object the name of your linear model object that you want to make a
#'summary table for.
#'@param titref the label name of the equation you made in Latex to cross
#'reference
#'@param labname the label name you want for this table
#'@param extracaption adds whatever text string you pass to the title of the
#'table.
#'@return xtable object
#'@seealso xtable
#'@export xtablelm
#'@examples
#'
#'##
#'
xtablelm <- function(lm.object, titref, labname, extracaption=NULL){
	require(xtable)
	
	x <- summary(lm.object)
	
	mat1 <-xtable(x, caption=paste("Summary Results for Regression in Equation \\eqref{",titref,"} ",extracaption,sep=""), label=labname, digits=c(0,3,3,2,4))
	
	# Computes and rounds off variables to be include at bottom of the table
	
	residse <- round(x$sigma,4)
	degf <- x$df[2]
	multr2 <- round(x$r.squared,5)
	adjr2 <- round(x$adj.r.squared,5)
	fstat <- round(x$fstatistic[1],3)
	fstatdf1 <- x$fstatistic[2]
	fstatdf2 <- x$fstatistic[3]
	fpval <- pf(x$fstatistic[1],x$fstatistic[2], x$fstatistic[3], lower.tail=FALSE)	
	# Adds the extra information to the end of the last row
	
	addtorow<- list()
	addtorow$pos <-list()
	addtorow$pos[[1]] <- c(dim(mat1)[1])
	addtorow$command[[1]]<-c(paste('\\hline \\multicolumn{5}{l}{Residual standard error:', residse, 'on', degf, 'degrees of freedom.} \\\\ \\multicolumn{5}{l}{Multiple R-squared:',multr2,'. Adjusted R-squared:', adjr2,'.} \\\\ \\multicolumn{5}{l}{F-statistic:', fstat, 'on', fstatdf1, 'and',  fstatdf2, 'DF, p-value:', signif(fpval, 4),'.} \\\\', sep=' '))
	
	
	print(mat1, add.to.row=addtorow, sanitize.text.function=NULL, caption.placement="top")

	#experiment with adding \hspace{1cm} to end of each row to improve spacing
	#for(i in 1:dim(mat1)[1]){
	#addtorow$pos[[i]] <- c(i)
	#}
	
	#for(i in 1:dim(mat1)[1]){
	#	if(i != dim(mat1)){addtorow$command[[i]]<-c(paste('\\hspace{1cm}'))}
	#	else{addtorow$command[[i]]<-c(paste('\\hline \\multicolumn{5}{l}{Residual standard error:', residse, 'on', degf, 'degrees of freedom.} \\\\ \\multicolumn{5}{l}{Multiple R-squared:', multr2,'.', 'Adjusted R-squared:', adjr2,'.','} \\\\ \\multicolumn{5}{l}{F-statistic:', fstat, 'on', fstatdf1, 'and',  fstatdf2, 'DF, p-value:', round(fpval, 6),'.','} \\\\', sep=' '))
	#	}
	#}
	#attr(mat1,"dimnames")[2]<-cat("Estimate","Std.Error","t value","Pr(>|t|)\\hspace{1cm}")

}



# Function will split column-wise a data.frame, matrix or a 2D array according
# to an INDEX (a factor) and apply a function if one supplied.
# 
# No simplification (sensu tapply) available at this time.
#
# Author: Roman Lustrik, May 2nd, 2011
###############################################################################




#'Split data over columns
#'
#'Split data column-wise on \code{data.frame}, \code{matrix} and \code{array}
#'or element-wise on a \code{list}.
#'
#'Function splits a \code{data.frame}, \code{matrix} and \code{array}
#'column-wise according to \code{INDEX} and \code{list} is sliced according to
#'\code{INDEX}. Output is returned as a list of the same length as the number
#'of levels in \code{INDEX}.
#'
#'@aliases splitc-package splitc
#'@param X A \code{data.frame}, \code{matrix}, \code{array} or a \code{list}.
#'@param INDEX A factor of \code{length(X)} (number of columns or list
#'elements). If not a factor, it will be coerced into one.
#'@param FUN A function to be applied to individual subset of data (each factor
#'level). If not provided (\code{NULL}), raw (split) data is returned.
#'@param \dots Additional arguments to \code{FUN}.
#'@return A list of the same length as there are factor levels in \code{INDEX}.
#'@note Simplification sensu \code{tapply} is not yet implemented.
#'@author Roman Lustrik \email{roman.lustrik@@biolitika.si}
#'@seealso \code{\link{tapply}}, \code{\link{by}}, \code{\link{aggregate}},
#'\code{\link{apply}}, \code{\link{split}}
#'@keywords manip
#'@export splitc
#'@examples
#'
#'my.list <- list(a = runif(5), b = runif(5), c = runif(5), d = runif(5), e = runif(10),
#'		f = runif(10), g = runif(10), h = runif(10), i = runif(10), j = runif(10))
#'my.df <- as.data.frame(my.list)
#'my.matrix <- as.matrix(my.df)
#'
#'ind <- factor(c(1,1,1,1, 2,3, 4,4,4,4))
#'ind2 <- factor(c(1,1,1,1, 2,3, 4,4,4,4), levels = 1:5)
#'
#'# Applies mean to each, you must use \code{colMeans}, 
#'#   as \code{mean} is deprecated for \code{data.frame}s
#'splitc(X = my.df, INDEX = ind, FUN = colMeans)
#'splitc(X = my.matrix, INDEX = ind2) # level 5 empty because not populated
#'splitc(X = my.list, INDEX = ind, FUN = sum) # applied to elements INDEX-wise 
#'
splitc <- function(X, INDEX, FUN = NULL, ...) {
  
  # Some initial checks
  wc <- class(X) # recognize class (JD, don't get funny ideas, wc is "what class?")
  
  if (!any(wc == c("matrix", "data.frame", "list", "array"))) 
    stop("Unrecognized class")
  
  FUN <- if (!is.null(FUN)) 
    match.fun(FUN)
  
  # Make some preparations before proceeding
  lvl <- as.list(levels(as.factor(INDEX))) # list of levels over which to iterate
  out.empty <- is.null(FUN) # if NULL, raw data should be returned
  
  
  # Extract columns for each level and if FUN supplied, apply a function
  # to the subset, otherwise return raw.
  out <- lapply(X = lvl, FUN = function(x, index, dta, fun, empty, ...) {
    
    my.subset <- index %in% x
    
    if (any(wc == "matrix" | wc == "data.frame" | wc == "array")) {
      pr <- dta[, my.subset, drop = FALSE]
      if (length(pr) < 1) return(NULL)
      else
        if (empty) return(pr) else return(do.call("FUN", list(pr, ...)))
    }
    
    if (wc == "list") {
      pr <- dta[my.subset]
      if (length(pr) < 1) return(NULL)
      else
        if (empty) return(pr) else return(lapply(X = pr, FUN = fun, ...))
    }
  }, index = INDEX, dta = X, fun = FUN, empty = out.empty, ...)
  
  return(out)
}
