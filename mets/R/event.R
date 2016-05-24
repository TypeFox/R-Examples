##' Constructur for Event History objects 
##'
##' ... content for details
##' @title Event history object
##' @param time Time 
##' @param time2 Time 2
##' @param cause Cause
##' @param cens.code Censoring code (default 0)
##' @param ... Additional arguments
##' @return Object of class Event (a matrix) 
##' @author Klaus K. Holst
##' @aliases Event rbind.Event
##' @examples
##' t1 <- 1:10
##' t2 <- t1+runif(10)
##' ca <- rbinom(10,2,0.4)
##' (x <- Event(t1,t2,ca))
Event <- function(time,time2=TRUE,cause=NULL,cens.code=0,...) {
    out <- cbind(time,time2,cause)
    if (!missing(cause)) {
        colnames(out) <- c("entry","exit","cause")
    } else {
        colnames(out) <- c("exit","cause")
    }
    class(out) <- "Event"
    attr(out,"cens.code") <- cens.code
    attr(out,"entry") <- ncol(out)==3
    return(out)
}

as.matrix.Event <- function(x,...) structure(x,class="matrix")

as.character.Event <- function(x,...) {
    if (ncol(x)==3) { 
        res <- paste("(",format(x[,1],...),";",
                     format(x[,2],...),":",
                     format(x[,3],...),"]",sep="")
    } else {
        res <- paste(format(x[,1],...),":",format(x[,2],...),sep="")
    }
    return(res)
}

format.Event <- function(x, ...) as.character.Event(x,...)

as.data.frame.Event <- as.data.frame.model.matrix

print.Event <- function(x,...) {
    print(as.matrix(x),...,quote=FALSE)
}

summary.Event <- function(object,...) {    
    print(table(object[,"cause"]))
    print(summary(object[,"exit"]))
}

"[.Event" <- function (x, i, j, drop = FALSE) 
{
    if (missing(j)) {
        atr <- attributes(x)
        class(x) <- "matrix"
        x <- x[i, , drop = FALSE]
        class(x) <- "Event"
        atr.keep <- c("cens.code","entry")
        attributes(x)[atr.keep] <- atr[atr.keep]
        x
    }
    else {
        class(x) <- "matrix"
        NextMethod("[")
    }
}

rbind.Event <- function(...) {
  dots <- list(...)
  cens.code <- attributes(dots[[1]])$cens.code
  ncol <- dim(dots[[1]])[2]
  nrow <- unlist(lapply(dots,nrow))
  cnrow <- c(0,cumsum(nrow))
  M <- matrix(ncol=ncol,nrow=sum(nrow))
  for (i in 1:length(dots)) {
    M[(cnrow[i]+1):cnrow[i+1],] <- dots[[i]]
  }
  x <- c(); for (i in 1:ncol(M)) x <- c(x,list(M[,i]))
  x <- c(x,list(cens.code=cens.code))
  do.call("Event",x)
} 


