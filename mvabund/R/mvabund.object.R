################################################################################
# function to create a mvabund object                                          #
################################################################################

mvabund <- function(..., row.names = NULL, check.rows = FALSE,
  check.names = TRUE, var.names=NULL, neg=FALSE, na.rm=FALSE){

  mc <- match.call(expand.dots= FALSE)
  mvabund.obj <- mc$...
  mvabund.obj <- lapply( mvabund.obj, eval, parent.frame() )
  
  if( length(mvabund.obj) ==1 & all(class(mvabund.obj[[1]])=="mvabund")) {
  class(mvabund.obj[[1]]) <- c("mvabund", "matrix")
  mvabund.obj <- data.frame(mvabund.obj[[1]], row.names=row.names, check.rows=check.rows,
      check.names=check.names)
  } else
    mvabund.obj <- data.frame(..., row.names=row.names, check.rows=check.rows,
      check.names=check.names)

names <- dimnames(mvabund.obj)
mvabund.obj <- as.matrix(mvabund.obj)
dimnames(mvabund.obj) <- names

  N <- nrow(mvabund.obj)

	if(na.rm & any(is.na(mvabund.obj))) {  
	mvabund.obj <- mvabund.obj[ -(which(is.na(mvabund.obj))%%N) , , drop=FALSE]
	N <- nrow(mvabund.obj)
	}
	
if (any(mvabund.obj[!is.na(mvabund.obj)]<0) & !neg )
     warning("There are some negative values in your abundance data, which is not expected if your data are counts - please check.")

 p <- dim(mvabund.obj)[2]

 if(!is.null(var.names)) {
  if(length(var.names)==p ) {
   if(check.names) {
	var.names <- make.names(var.names, unique = TRUE)  }
  	dimnames(mvabund.obj)[[2]] <- var.names
  } else
  stop("length of 'row.names' not equal to number of variables")
 }

class(mvabund.obj) <- c("mvabund","matrix")
if(!is.numeric(mvabund.obj))
    warning("The object is not numeric. Further analysis will fail.")

return(mvabund.obj)
}



############ a function to turn its argument into a mvabund object #############
as.mvabund <- function(x) {
if (is.mvabund(x)) x else mvabund(x)
}


########### a function to test if its argument is a mvabund object #############
is.mvabund <- function(x) {
inherits(x, "mvabund")
}



# setMethod("show", "mvabund",
# function(object) {
# uncl <- unclass(object)
# show(uncl)
# })   


print.mvabund <- function(x,...) {
uncl <- unclass(x)
print(uncl,...)
}

# setMethod("print", "mvabund",
# function(x,...) {
# uncl <- unclass(x)
# print(uncl,...)
# })


