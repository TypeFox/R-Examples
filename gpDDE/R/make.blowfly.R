
make.blowfly <- function(){
    fn <- function(t, y, p, more){
        r = y
        y.d <- more$y.d[,1]
        r[,"y"] <- p["c"] * y.d * exp(-y.d / p["N0"]) - p["a"] * y
        return(r)
    }

    dfdx <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        r[,"y", "y"] <- -p["a"]
        return(r)
    }

    dfdx.d <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        y.d <- more$y.d[,1]
        r[,"y", "y"] <- p["c"] * (1 - y.d / p["N0"]) * exp(-y.d / p["N0"])
        return(r)
    }

    dfdp <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), names(p))
        y.d <- more$y.d[,1]
        r[,"y", "c"] <- y.d * exp(-y.d / p["N0"])
        r[, "y", "N0"] <- p["c"] * (y.d / p["N0"])^2 * exp(-y.d / p["N0"])
        r[, "y", "a"] <- -y
        return(r)
    }

    d2fdx2 <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    d2fdxdp <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
        y.d <- more$y.d[,1]
        r[,"y", "y", "a"] <- -1
        return(r)
    }

    d2fdx.ddp <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
        y.d <- more$y.d[,1]
        r[,"y", "y", "c"] <- exp(-y.d / p["N0"]) * (1 - y.d / p["N0"])
        r[,"y", "y", "N0"] <- p["c"] * exp(-y.d / p["N0"]) * (2 * y.d / p["N0"]^2 - y.d^2 / p["N0"]^3)
        return(r)
    }

    d2fdx.d2 <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        y.d <- more$y.d[,1]
        r[,"y","y","y"] <- p["c"] * exp(-y.d / p["N0"]) * (y.d / p["N0"]^2 - 2 / p["N0"])
        return(r)
    }

    d2fdxdx.d <- function(t, y, p, more){
        r <- array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
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


