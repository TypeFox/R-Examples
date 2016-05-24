#' @title Make Vector Disease functions
#'
#' @return A list of functions that calculate the derivatives of vector disease model.
#' @export
#' @examples
#' make.vector.disease.fn()
#' vector.disease.fn <- make.vector.disease.fn()
make.vector.disease.fn <- function(){
    fn <- function (times, y, p, more){
        y.d <- more$y.d[,1]
        r <- y
        dimnames(r) <- dimnames(y)
        r[, "S"] <- (p["b"] + sin(times))* y.d * (1 - y[,"S"]) - p["a"] * y[,"S"]
        return(r)
    }

    dfdx <- function(times, y, p, more){
        y.d <- more$y.d[,1]
        r = array(0, c(dim(y), dim(y)[2]))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        r[,"S", "S"] <- -(p["b"] + sin(times)) * y.d - p["a"]
        return(r)
    }

    dfdp <- function(times, y, p, more){
        y.d <- more$y.d[,1]
        r = array(0, c(dim(y), 2))
        dimnames(r) = list(NULL, colnames(y), names(p)[c(1,2)])
        r[, "S", "a"] <- - y[,"S"]
        r[, "S", "b"] <- y.d * (1 - y[,"S"])
        return(r)
    }

    d2fdxdp <- function (times, y, p, more){
        y.d <- more$y.d[,1]
        r = array(0, c(dim(y), dim(y)[2], 2))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p)[c(1,2)])
        r[,"S", "S", "a"] <- -1
        r[,"S", "S", "b"] <- - y.d
        return(r)
    }

    d2fdx2 <- function(times, y, p, more){
        y.d <- more$y.d[,1]
        r = array(0, c(dim(y), 1, 1))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    d2fdx.ddp <- function(times, y, p, more){
        y.d <- more$y.d[,1]
        r <- array(0, c(dim(y), 1, 2))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p)[c(1:2)])
        r[,"S", "S", "a"] <- 0
        r[,"S", "S", "b"] <- 1 - y
        return(r)
    }


    dfdx.d <- function(times, y, p, more){
        y.d <- more$y.d[,1]
        r = array(0, c(dim(y), dim(y)[2]))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        r[,"S", "S"] <- (p["b"] + sin(times)) * (1 - y[,"S"])
        return(r)
    }

    d2fdxdx.d <- function(times, y, p, more){
        y.d <- more$y.d[,1]
        r = array(-(p["b"] + sin(times)), c(dim(y), 1, 1))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    d2fdx.d2 <- function(times, y, p, more){
        y.d <- more$y.d[,1]
        r = array(0, c(dim(y), 1, 1))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    return(list(
        fn = fn, dfdx = dfdx,
        dfdp = dfdp, d2fdx2 = d2fdx2,
        d2fdxdp = d2fdxdp, d2fdx.ddp = d2fdx.ddp,
        dfdx.d = dfdx.d, d2fdx.ddx = d2fdxdx.d,
        d2fdxdx.d = d2fdxdx.d, d2fdx.d2 = d2fdx.d2
        ))
}


