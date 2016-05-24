###############################################################################
## Project: ETLUtils
## Content: Fast merging based on match
##
## Author: jan
## Creation date: Mar 25, 2012, 2:56:02 PM
## File: matchmerge.R
###############################################################################


#' Merge two data frames (fast) by common columns by performing a left (outer) join or an inner join.  
#' 
#' Merge two data frames (fast) by common columns by performing a left (outer) join or an inner join.\cr
#' The data frames are merged on the columns given by by.x and by.y. Columns can be specified only by name.
#' This differs from the merge function from the base package in that merging is done based on 1 column key only. 
#' If more than one column is supplied in by.x and by.y, these columns will be concatenated together 
#' to form 1 key which will be used to match. 
#' Alternatively, by.x and by.y can be 2 vectors of length NROW(x) which will be used as keys.\cr
#' 
#' The rows in the right hand side data frame that match on the specific key are extracted, and joined together
#' with the left hand side data frame.\cr
#' 
#' Merging is done based on the match function on the key value.
#' This makes the function a lot faster when compared to applying merge, especially for large data frames (see the example).
#' And also the memory consumption is a lot smaller.\cr
#' 
#' In SQL database terminology, the default value of all.x = FALSE gives a natural join, 
#' a special case of an inner join. Specifying all.x = FALSE gives a left (outer) join. 
#' Right (outer) join or (full) outer join are not provided in this function. 
#' 
#' @param x the left hand side data frame to merge
#' @param y the right hand side data frame to merge \cr
#' or a vector in which case you always need to supply by.y as a vector, make sure
#' by.iskey is set to TRUE and provide in add.columns the column name for which y will be relabelled to in the joined data frame (see the example). 
#' @param by.x either the name of 1 column in x or a character vector of length NROW(x) 
#' which will be used as key to merge the 2 data frames
#' @param by.y either the name of 1 column in y or a character vector of length NROW(x) 
#' which will be used as key to merge the 2 data frames. 
#' Duplicate values in by.y are not allowed.
#' @param all.x logical, if TRUE, then extra rows will be added to the output, one for each row in x that has no matching row in y. 
#' These rows will have NAs in those columns that are usually filled with values from y. 
#' The default is FALSE, so that only rows with data from both x and y are included in the output. The
#' default value corresponds to an inner join. If TRUE is supplied, this corresponds to a left (outer) join.
#' @param by.iskey Logical, indicating that the by.x and the by.y inputs are vectors of length NROW(x) and NROW(y)
#' instead of column names in x and y. If this is FALSE, the input columns will be pasted together to create a key
#' to merge upon. Otherwise, the function will use the by.x and by.y vectors directly as matching key. 
#' Defaults to FALSE indicating the by.x and by.y are column names in x and y.
#' @param suffix a character string to be used for duplicate column names in x and y to make the y columns unique.
#' @param add.columns character vector of column names in y to merge to the x data frame. Defaults to all columns in y.
#' @param check.duplicates checks if by.y contains duplicates which is not allowed. Defaults to TRUE.
#' @param trace logical, indicating to print some informative messages about the progress
#' @return data frame with x joined with y based on the supplied columns.
#' The output columns are the columns in x followed by the extra columns in y. 
#' @export
#' @seealso \code{\link{cbind}, \link{match}, \link{merge}}
#' @examples
#' left <- data.frame(idlhs = c(1:4, 3:5), a = LETTERS[1:7], stringsAsFactors = FALSE)
#' right <- data.frame(idrhs = c(1:4), b = LETTERS[8:11], stringsAsFactors = FALSE)
#' ## Inner join
#' matchmerge(x=left, y=right, by.x = "idlhs", by.y = "idrhs")
#' 
#' ## Left outer join in 2 ways
#' matchmerge(x=left, y=right, by.x = "idlhs", by.y = "idrhs", all.x=TRUE)
#' matchmerge(x=left, y=right, by.x = left$idlhs, by.y = right$idrhs, all.x=TRUE, by.iskey=TRUE)
#' 
#' ## Show usage when y is just a vector instead of a data.frame
#' matchmerge(x=left, y=right$b, by.x = left$idlhs, by.y = right$idrhs, all.x=TRUE, 
#' by.iskey=TRUE, add.columns="b.renamed")
#' 
#' ## Show speedup difference with merge
#' \dontrun{
#' size <- 100000 
#' dimension <- seq(Sys.Date(), Sys.Date()+10, by = "day")
#' left <- data.frame(date = rep(dimension, size), sales = rnorm(size))
#' right <- data.frame(date = dimension, feature = dimension-7, feature = dimension-14)
#' dim(left)
#' dim(right)
#' print(system.time(merge(left, right, by.x="date", by.y="date", all.x=TRUE, all.y=FALSE)))
#' print(system.time(matchmerge(left, right, by.x="date", by.y="date", all.x=TRUE, by.iskey=FALSE)))
#' }
#' ## Show example usage 
#' products <- expand.grid(product = c("Pepsi", "Coca Cola"), type = c("Can","Bottle"), 
#' size = c("6Ml","8Ml"), distributor = c("Distri X","Distri Y"), salesperson = c("Mr X","Mr Y"), 
#' stringsAsFactors=FALSE)
#' products <- products[!duplicated(products[, c("product","type","size")]), ]
#' products$key <- paste(products$product, products$type, products$size, sep=".")
#' sales <- expand.grid(item = unique(products$key), sales = rnorm(10000, mean = 100))
#' str(products)
#' str(sales)
#' info <- matchmerge(x=sales, y=products, 
#'   by.x=sales$item, by.y=products$key, all.x=TRUE, by.iskey=TRUE, 
#'   add.columns=c("size","distributor"), check.duplicates=FALSE)
#' str(info)
#' tapply(info$sales, info$distributor, FUN=sum)
matchmerge <- function(x, y, by.x, by.y, all.x=FALSE, by.iskey=FALSE, suffix = ".y", add.columns=colnames(y), check.duplicates=TRUE, trace=FALSE){
	if(!"data.frame" %in% class(x)){
		stop("x should be of class data.frame")
	}
	if(!"data.frame" %in% class(y)){
		if(!is.vector(y) | by.iskey == FALSE){
			stop("y should be of class data.frame or a vector in which case by.iskey should be TRUE and you should supply by.y as a vector")	
		}else{
			if(length(y) != length(by.y)){
				stop("y is a vector so you have to make sure y and by.y are of the same length")
			}
		}
	}else{
		if(sum(!add.columns %in% colnames(y))){
			stop("all add.columns should be in colnames(y)")
		}	
	}
	
	## Make the keys
	if(by.iskey == FALSE){
		if(trace) message(sprintf("%s: Creating key variables", Sys.time()))
		if(length(by.x) > 1){
			by.lhs.vector <- do.call(paste, as.list(x[, by.x]))	
		}else{
			by.lhs.vector <- x[[by.x]]
		}
		if(length(by.y) > 1){
			by.rhs.vector <- do.call(paste, as.list(y[, by.y]))	
		}else{
			by.rhs.vector <- y[[by.y]]
		}
		#add.columns <- add.columns[!add.columns %in% by.x]
	}else{
		if(length(by.x) != NROW(x)){
			stop("by.x is not of length NROW(x)")
		}
		if(length(by.y) != NROW(y)){
			stop("by.y is not of length NROW(y)")
		}
	}	
	## Remove data from x if inner join
	if(all.x == FALSE){
		if(trace) message(sprintf("%s: Removing left hand side data for inner join", Sys.time()))
		if(by.iskey == FALSE){
			idx <- which(by.lhs.vector %in% by.rhs.vector)
			by.lhs.vector <- by.lhs.vector[idx]
		}else{
			idx <- which(by.x %in% by.y)
		}	
		x <- x[idx, , drop=FALSE]
	}
	## Find the overlapping keys 
	if(trace) message(sprintf("%s: Searching for key matches", Sys.time()))
	if(by.iskey == FALSE){
		if(check.duplicates == TRUE){
			if(sum(duplicated(by.rhs.vector)) > 0){
				stop("Key of y contains doubles which is not allowed for this inner or left outer join")
			}	
		}
		idx <- match(x=by.lhs.vector, table = by.rhs.vector)
	}else{		
		if(check.duplicates == TRUE){
			if(sum(duplicated(by.y)) > 0){
				stop("Key of y contains doubles which is not allowed for this inner or left outer join")
			}	
		}
		idx <- match(x=by.x, table = by.y)
	}
	## Join the data
	if(trace) message(sprintf("%s: Joining the data frames", Sys.time()))
	if("data.frame" %in% class(y)){
		overlapping.columns <- which(add.columns %in% colnames(x))
		if(length(overlapping.columns) > 0){
			result <- cbind(x, renameColumns(
							y[idx, add.columns, drop=FALSE], 
							from = add.columns[overlapping.columns], 
							to = paste(add.columns[overlapping.columns], suffix, sep="")))
		}else{
			result <- cbind(x, y[idx, add.columns, drop=FALSE])
		}
	}else{
		result <- cbind(x, y[idx])
		colnames(result)[ncol(result)] <- add.columns[1]
	}
	
	rownames(result) <- rownames(x)
	result
}


