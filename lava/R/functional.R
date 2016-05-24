##' @export
"functional<-" <- function(x,...,value) UseMethod("functional<-")

##' @export
"functional<-.lvm" <- function(x,to,from,...,value) {
    if (inherits(to,"formula")) {
        yy <- decomp.specials(getoutcome(to))
        ##xx <- attributes(terms(to))$term.labels
        myvars <- all.vars(to)
        xx <- setdiff(myvars,yy)
        if (length(yy)*length(xx)>length(value) & length(value)!=1) stop("Wrong number of values")
        count <- 0
        for (y in yy) {
            count <- count+1
            for (i in seq_along(xx)) {
                suppressWarnings(x <- regression(x,to=y,from=xx[i],silent=TRUE))
                count <- count+1
                if (length(value)==1) {
                    functional(x, to=y, from=xx[i],...) <- value
                } else
                    functional(x, to=y, from=xx[i],...) <- value[[count]]
            }
        }
        return(x)
    }

    if (missing(from) | missing(to))
        return(x)

    edges <- paste(from,to,sep="~")
    x$attributes$functional[[edges]] <- value
    return(x)
}

##' @export
"functional" <- function(x,...) UseMethod("functional")

##' @export
functional.lvm <- function(x,to,from,...) {
    if (missing(from))
        return(x$attributes$functional)

    edges <- paste(from,to,sep="~")
    x$attributes$functional[edges]
}
