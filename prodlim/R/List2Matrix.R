### List2Matrix.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 21 2015 (07:01) 
## Version: 
## last-updated: Sep 29 2015 (06:32) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: Reduce a list to a matrix or data.frame and add list names as new columns
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' This function is used by summary.prodlim to deal with results.
##'
##' Reduction is done with rbind.
##' @title Reduce list to a matrix or data.frame with names as new columns
##' @param list A named list which contains nested lists 
##' @param depth The depth in the list hierarchy until an rbindable object
##' @param names Names for the list variables
##' @return Matrix or data.frame.
##' @examples
##' 
##' x=list(a=data.frame(u=1,b=2,c=3),b=data.frame(u=3,b=4,c=6))
##' List2Matrix(x,depth=1,"X")
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
List2Matrix <- function(list,depth,names){
    if (missing(names)) names <- paste0("D",1:depth)
    switch(as.character(depth),
           "1"={
               dims <- lapply(list,dim)
               cols <- sapply(dims,function(x)x[[2]])
               rows <- sapply(dims,function(x)x[[1]])
               stopifnot(length(unique(cols))==1)
               nl <- names(list)
               M <- do.call("rbind",list)
               rownames(M) <- NULL
               M <- cbind(rep(nl,rows),M)
               colnames(M)[1] <- names[1]
               M},
           "2"={
               List2Matrix(lapply(list,List2Matrix,depth=1,names=names[2]),
                           depth=1,
                           names=names[1])},
           "3"={
               List2Matrix(lapply(list,function(l){
                                      List2Matrix(lapply(l,List2Matrix,depth=1,names[3]),
                                                  depth=1,
                                                  names=names[2])
                                  }), depth=1,names=names[1])},
           stop("Cannot do this depth."))
}
#----------------------------------------------------------------------
### List2Matrix.R ends here
