## A helper function to reorder vector v (if named) into order
## specified by names. Copied from ergm::ergm.model.utils.R.
vector.namesmatch<-function(v,names,errname=NULL){
  if(is.null(errname)) errname <- deparse(substitute(v))

  if (is.null(names(v))){
    if(length(v) == length(names)){
      names(v) <- names
    }else stop('Length of "', errname, '" is ', length(v), " but should be ", length(names),".")
  }else{
    if(length(v) == length(names)
       && length(unique(names(v)))==length(v)
       && length(unique(names))==length(names)
       && all(sort(names(v)) == sort(names))){
      namesmatch <- match(names, names(v))
      v <- v[namesmatch]
    }else stop('Name missmatch in "', errname,'". Specify by position.')
  }
  v
}

## Compress a data frame by eliminating duplicate rows while keeping
## track of their frequency.
compress.data.frame<-function(x){
  x<-sort(x)
  firsts<-which(!duplicated(x))
  freqs<-diff(c(firsts,nrow(x)+1))
  x<-x[firsts,]
  list(rows=x,frequencies=freqs)
}

## Sorts rows of a data frame in lexicographic order.
sort.data.frame<-function(x, decreasing=FALSE, ...){
  x[do.call(order,c(sapply(seq_along(x),function(i)x[[i]],simplify=FALSE), decreasing=decreasing)),]
}

## Return the first non-NULL argument. If all arguments are NULL, return NULL.
NVL <- function(...){
  for(x in list(...))
    if(!is.null(x)) break
  x
}

## Return the first non-try-error argument. If all arguments are try-errors, stop with an error.
ERRVL <- function(...){
  for(x in list(...))
    if(!inherits(x, "try-error")) return(x)
  stop("No non-error expressions passed.")
}

## Only run expr if environment variable testvar is set to specified values. Otherwise, skip them and optionally print a message documenting this.
opttest <- function(expr, testname=NULL, testvar="ENABLE_statnet_TESTS", yesvals=c("y","yes","t","true","1"), lowercase=TRUE){
  testval <- Sys.getenv(testvar)
  if(lowercase) testval <- tolower(testval)
  if(testval %in% yesvals)
    eval.parent(expr)
  else
    if(!is.null(testname))
      message(testname," test(s) skipped. Set ",testvar," environment variable to run.")
}
