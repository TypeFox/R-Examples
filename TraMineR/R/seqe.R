## ========================================
## Check if subsequence
## ========================================
is.seqe<-function(s){
#   return(.Call("istmrsequence",s, PACKAGE="TraMineR"))
    return(inherits(s,"seqe"))
}
is.seqelist<-function(s){
#   return(.Call("istmrsequence",s, PACKAGE="TraMineR"))
   return(inherits(s,"seqelist"))
}
#as.seqelist<-function(s){
#   return(.Call("istmrsequence",s, PACKAGE="TraMineR"))
#    class(s)<-c("seqelist","list")
#    return(s)
#}

###Methods taken from survival package

"[.seqelist" <- function(x, i,j, drop=FALSE) {
    # If only 1 subscript is given, the result will still be a Surv object
    #  If the second is given extract the relevant columns as a matrix
  if (missing(j)) {
  	temp <- class(x)
  	type <- attr(x, "type")
  	class(x) <- NULL
  	x <- x[i, drop=FALSE]
  	class(x) <- temp
  	attr(x, "type") <- type
  	x
	} else {
  	class(x) <- NULL
  	NextMethod("[")
	}
}
Math.seqelist <- function(...){
  stop("Invalid operation on event sequences")
}
Ops.seqelist  <- function(...){
  stop("Invalid operation on event sequences")
}
Summary.seqelist<-function(...) {
  stop("Invalid operation on event sequences")
}
Math.seqe <- function(...)  {
  stop("Invalid operation on event sequences")
}
Ops.seqe  <- function(...)  {
  stop("Invalid operation on event sequences")
}
Summary.seqe<-function(...) {
  stop("Invalid operation on event sequences")
}

levels.seqe<-function(x,...){
  if(!is.seqe(x))stop("x should be a seqe object. See help on seqecreate.")
  return(.Call(TMR_tmrsequencegetdictionary,x))
}

levels.seqelist<-function(x,...){
  if(!is.seqelist(x))stop("x should be a seqelist. See help on seqecreate.")
  if(length(x)>0) return(.Call(TMR_tmrsequencegetdictionary,x[[1]]))
}
## ========================================
## Return a string representation of a sequence
## ========================================

str.seqelist<-function(object,...){
#message("Event sequence analysis module is still experimental")
  if(is.seqelist(object)){
      object<-cat(as.character(object),"\n")
  }else if (is.seqe(object)){
    object<-cat(as.character(object),"\n")
  }else{
    stop("object should be a seqelist. See help on seqecreate.")
  }
  NextMethod("str")
}
str.seqe<-function(object,...){
#  seqestr(s)
  if(!is.seqe(object))stop("object should be a seqe object. See help on seqecreate.")
  object <- .Call(TMR_tmrsequencestring, object)
  NextMethod("str")
}

as.character.seqe<-function(x, ...){
#  seqestr(s)
  if(!is.seqe(x))stop("x should be a seqe object. See help on seqecreate.")
  x<-.Call(TMR_tmrsequencestring,x)
  NextMethod("as.character")
}

as.character.seqelist<-function(x, ...){
  tmrsequencestring.internal<-function(s){
    if(is.seqe(s)){
      return(.Call(TMR_tmrsequencestring, s))
    }
    return(as.character(s))
  }
  if(!is.seqelist(x))  stop("x should be a seqelist object. See help on seqecreate.")
  x <- as.character(sapply(unlist(x), tmrsequencestring.internal))
  NextMethod("as.character")
}

## ========================================
## Print sequences
## ========================================

#seqeprint<-function(s){
#  print(seqestr(s))
#}
print.seqe<-function(x,quote = FALSE, ...){
  x <- as.character(x)
  print(x, quote=quote, ...)
}
print.seqelist<-function(x, quote = FALSE, ...){
  x <- as.character(x)
  print(x, quote=quote, ...)
}

## ========================================
## Plot sequences
## ========================================

plot.seqe <- function(x, type = "pc", ...) {
  if (type == "pc") {
    seqpcplot(x, ...)
  }
}

plot.seqelist <- function(x, type = "pc", ...) {
  if (type == "pc") {
    seqpcplot(x, ...)
  }
}
