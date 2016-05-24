###############################################################################
# [.mvabund: a function to use [ on mvabund                                   #
# this function should be invisible, should not be used directly,             #
# only through [                                                              #
###############################################################################

'[.mvabund' <-  function(x,i, j=NULL, ... , drop = TRUE ) {
if(!is.mvabund(x)) stop("use only with a mvabund object")

if(is.null(dim(x)) | length(dim(x))>2) x <- mvabund(c(x))

dots <- match.call(expand.dots = FALSE)$...  

classx <- class(x)

if (length(class(x))>1){
# x is a mvabund and a matrix object.

    n <- NROW(x)
    p <- NCOL(x)
    x <- matrix(x, nrow=n, ncol=p, dimnames = dimnames(x))
    
    if (missing(i)) i <- 1:n
    
    if(p>1 & length(i)==n*p)  {
    # If only one index argument is given, for all the values of x.
    
        if(length(dots)>0) {
          b=x[i,..., drop=drop] } else
          b=x[i, drop=drop] 
    
    } else { 
    
    if (missing(j)) j <- 1:p
    if(length(dots)>0) {
        b=x[i,j,..., drop=drop] } else
        b=x[i,j, drop=drop] 
    
    }
 
} else if (length(class(x))==1){
# x is a mvabund object.

    n <- NROW(x)
    p <- NCOL(x)
    x <- matrix(x, nrow=n, ncol=p, dimnames = dimnames(x) )
    if (missing(i)) i <- 1:n
    
    if(length(i)==n*p)  {
    # If only one index argument is given, for all the values of x.

        if(length(dots)>0) {
          b=x[i,..., drop=drop] } else
          b=x[i, drop=drop] 
    
    } else {    
        if (missing(j)) j <- 1:p
        if(length(dots)>0)
            b=x[i,j,..., drop=drop] else
            b=x[i,j, drop=drop]
    }     
} else stop("class is missing")
classb <- class(b) 

if (is.matrix(b) | length(dim(array(b)))>=2) {
  # If b is a matrix or an array it should additonally be a mvabund object.
    class(b) <- c("mvabund", c(classb))
} else {
  # The extracted object should only be mvabund, when it is multidimensional.
  class(b) <- classb[classb != "mvabund"]
}

return(b)

}

