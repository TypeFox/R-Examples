##' Integrate characters with careful ordering of result
##'
##'
##' @param x character vector
##' @param y character vector
##' @importFrom stringr str_split_fixed
##' @return union of x and y, with
##' \emph{write something about expected return ordering}.
vectorIntegrate <- function(x, y) {
    vectorIntegrate.eval(x, y, vectorIntegrate.indices(x, y))
}

##' Test equality and order of a vector
##'
##'
##' @param x character
##' @param y character
##' @param order logical. stricter test for order of x and y
##' @return logical
vecEqual <- function (x, y, order=FALSE) {
	if (length(x) != length(y)) {
		eq <- FALSE
	} else if (order) {
		eq <- TRUE
		i <- 0
		while (eq & i<length(x)) {
			i <- i+1
			eq <- eq & x[i]==y[i]
		}
	} else {
		eq <- setequal(x,y)
	}
	return(eq)
}

rowIntegrate <- function (x, y, sep="QQQQQ") {
    stopifnot(ncol(x)==ncol(y))
    n <- ncol(x)
    x <- apply(x, 1, paste, collapse=sep)
    y <- apply(y, 1, paste, collapse=sep)

    out <- vectorIntegrate(x, y)
    return(str_split_fixed(out, sep, n))
}



vectorIntegrate.indices <- function (x, y) {
    if (length(x)==0) {
        if (length(y)==0) return(data.frame(vector=numeric(0), ind=numeric(0)))
        return(data.frame(vector=2, ind=1:length(y)))
    }

    if (length(y)==0) {
        return(data.frame(vector=1, ind=1:length(x)))
    }

    if (vecEqual(x,y)) {
        return(data.frame(vector=3, ind=1:length(x)))
    }

    if (length(intersect(x, y))==0) {
        return(rbind(data.frame(vector=1, ind=1:length(x)),
            data.frame(vector=2, ind=1:length(y))))
    }

    if (length(x)>length(y)) {
        a <- x
        x <- y
        y <- a
        flipped <- TRUE
    } else {
        flipped <- FALSE
    }

    out <- data.frame(vector=3, ind=match(y, x))
    if (any(is.na(out$ind))) {
        out[is.na(out$ind),] <- data.frame(vector=2, ind=which(is.na(out$ind)))
    }

    adds <- !(1:length(x) %in% out$ind[out$vector==3])
    #adds <- x %in% t1

    while (any(adds)) {
        these <- firstrun(adds)
        this <- which(out$vector %in% c(1,3) && out$ind %in% (max(these)+1))
        if (length(this)==0) {
            out <- rbind(out, data.frame(vector=1, ind=these))
        } else {
            t2 <- rbind(data.frame(vector=1, ind=these), out[this,])
            if (this>1) t2 <- rbind(out[1:(this-1),], t2)
            if (this<nrow(out)) t2 <- rbind(t2, out[(this+1):nrow(out),])
            out <- t2
        }

        adds <- !(1:length(x) %in% out$ind[out$vector %in% c(1,3)])
    }

    if (flipped) {
        out$vector[out$vector<3] <- 3-out$vector[out$vector<3]
        if (any(out$vector==3)) {
            out$ind[out$vector==3] <- which(y %in% x[out$ind[out$vector==3]])
        }
    }
    return(out)

}

vectorIntegrate.eval <- function (x, y, indices, joinFunction=NULL) {
    z <- rep(NA, nrow(indices))
    if (length(z)) {
        z[indices$vector %in% c(1,3)] <- x[indices$ind[indices$vector %in% c(1,3)]]
        z[indices$vector==2] <- y[indices$ind[indices$vector==2]]
    }
    return(z)
}

## returns the indices of the first run of TRUEs in a vector of logicals
firstrun <- function (x) {
    if (!any(x)) return(c())
    x <- which(x)
    out <- x[1]
    counter <- TRUE
    i <- out+1
    while (counter) {
        if (i %in% x) {
            out <- c(out, i)
            i <- i+1
        } else {
            counter <- FALSE
        }
    }
    return (out)
}

##' Integrate matrices based on their dim names
##'
##'
##' @param x matrix
##' @param y matrix
##' @param x.names a list of length 2
##' @param y.names a list of length 2
##' @param joinFunction function used to join things
##' @return union of x and y, with
##' \emph{write something about expected return ordering}.
matrixIntegrate <- function (x, y, x.names=dimnames(x), y.names=dimnames(y),
                            joinFunction=pasteExp) {

    ## if fn is NULL, just return the first one
    if (is.null(joinFunction)) joinFunction <- function (...) list(...)[[1]]

    row.ind <- vectorIntegrate.indices(x.names[[1]], y.names[[1]])
    col.ind <- vectorIntegrate.indices(x.names[[2]], y.names[[2]])
    r <- do.call("cbind", sapply(1:nrow(col.ind), function (z) {
        ## by column, integrate using row.ind
        this <- col.ind$vector[z]
        this.r <- rep(NA, nrow(row.ind))
        if (this==1) {
            this.r[row.ind$vector %in% c(1,3)] <- x[row.ind$ind[row.ind$vector %in% c(1,3)],
                col.ind$ind[z]]
        } else if (this==2) {
            this.r[row.ind$vector==2] <- y[row.ind$ind[row.ind$vector==2],
                col.ind$ind[z]]
            # ## need to know x's ind if 3:
            zbis <- match(x.names[[1]][row.ind$ind[row.ind$vector==3]],
                y.names[[1]])
            this.r[row.ind$vector==3] <- y[zbis,col.ind$ind[z]]
        } else if (this==3) {
            ## need to know y's col ind if 3:
            zbis <- which(y.names[[2]] %in% x.names[[2]][col.ind$ind[z]])
            this.r <- vectorIntegrate.eval(x[,col.ind$ind[z]], y[,zbis],
                row.ind)
            ## need to know y's row ind if 3:
            zrow <- match(x.names[[1]][row.ind$ind[row.ind$vector==3]],
                y.names[[1]])
            this.r[row.ind$vector==3] <- joinFunction(this.r[row.ind$vector==3],
                y[zrow,zbis])
        }
        return(this.r)
    }, simplify=FALSE))
    r.names <- list(vectorIntegrate.eval(x.names[[1]], y.names[[1]], row.ind),
                    vectorIntegrate.eval(x.names[[2]], y.names[[2]], col.ind))

    return(list(object=r, dimnames=r.names))
}

### This could go wrong in a number of places.
pasteExp <- function (..., sep="  |  ") {
    dots <- list(...)
    if (any(unlist(lapply(dots, length))>1)) {
        return(mapply("pasteExp", ..., USE.NAMES=FALSE))
    } else {
        dots <- unique(unlist(lapply(list(...), function (x) {
                x <- unlist(strsplit(x, sep, fixed=TRUE))
                if (length(x)>1) x <- sub("^\\(", "", sub("\\)$", "", x))
                return(x)
            })))
        if (length(dots)>1) {
            dots <- sapply(unique(dots),
                function (x) paste("(", x, ")", sep=""))
        }
        return(paste(dots, collapse=sep))
    }
}
