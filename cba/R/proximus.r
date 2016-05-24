
#
#   proximus.r - according to the paper:
#
#   M. Koyutürk, A. Graham, and N. Ramakrishnan. Compression, 
#   Clustering, and Pattern Discovery in Very High-Dimensional 
#   Descrete-Attribute Data Sets. IEEE Transactions On Knowledge 
#   and Data Engineering, Vol. 17, No. 4, (April) 2005
#
#   Contents: 
#
#   wrapper(s) for my C implementation of PROXIMUS. a better 
#   implementation may use two sparse matrices holding the pair 
#   of approximating matrices X and Y so that hat A = X * Y.
#
#   Version: 0.1-1
#
#   (C) ceeboo, 2005

proximus <- function(x, max.radius=2, min.size=1, min.retry=10, max.iter=16, 
                                                                debug=FALSE) {
    if (!is.logical(x))
       stop(paste(sQuote("x"),"not logical"))

    storage.mode(max.radius) <- storage.mode(min.size) <- "integer"
    storage.mode(min.retry) <- storage.mode(max.iter) <- "integer"
    storage.mode(debug) <- "logical"
    
    obj <- .Call(R_proximus, x, max.radius, min.size, min.retry, max.iter, 
                                                                 debug)
    obj$max.radius <- max.radius
    obj$min.size <- min.size
    obj$rownames <- rownames(x)
    obj$colnames <- colnames(x)
    class(obj) <- c("proximus")
    invisible(obj)
}

# get the full storage representation + pattern (cluster) labels

fitted.proximus <- function(object, drop=FALSE, ...) {
    x <- matrix(FALSE, nrow=object$nr, ncol=object$nc)
    c <- vector("integer", object$nr)
    for (i in 1:length(object$a)) {
        x[object$a[[i]]$x, object$a[[i]]$y] <- TRUE
        c[object$a[[i]]$x] <- i
    }
    k <- rep(TRUE, object$nr)   # keep
    if (drop) {
       for (i in 1:length(object$a)) 
           if (length(object$a[[i]]$x) < object$min.size ||
               object$a[[i]]$r > object$max.radius)
              k[object$a[[i]]$x] <- FALSE
       x <- x[k,]
       c <- c[k]
    }
    rownames(x) <- object$rownames[k]
    colnames(x) <- object$colnames
    attr(c, "Index") <- which(k)
    #
    x <- list(x=x, pl=factor(c))
    x
}

###

print.proximus <- function(x, ...) {
    cat("an object of class:",class(x),"\n")
    invisible(x)
}

summary.proximus <- function(object, ...) {
    n <- length(object$a)
    s <- as.data.frame(matrix(nrow=n, ncol=7))
    names(s) <- c("Size","Length","Radius","Error","Fnorm","Jsim","Valid")
    e <- j <- 0
    for (i in 1:n) {                                # pattern summaries
        a  <- object$a[[i]]                         # approximation
        nx <- length(a$x)
        ny <- length(a$y)
        s[i,] <- c(nx, ny, a$r, 
                   (a$n - a$c) / (nx * object$nc),  # Error
		           sqrt(a$n - a$c),                 # Frobenius norm
                   if (a$c == 0 && ny == 0) 
                      1                             # definition!
                   else
                      1 / (1 + 2 * (a$n - a$c) /
                                   (a$c + nx * ny)),# Jaccard
                   (nx  >= object$min.size &
                    a$r <= object$max.radius))      # valid
        e <- e + a$n - a$c			                # total Error
        j <- j + a$c + nx * ny			            # total Jaccard
    }
    storage.mode(s[,7]) <- "logical"
    s <- list(nr=object$nr, 
              nc=object$nc,
              error=e / (object$nr * object$nc), 
              fnorm=sqrt(e),
              jsim=if (j == 0 && e == 0) 
                      1                             # definition!
                   else                  
                      j / (j + e / 2),
              valid=sum(s$Valid),
              pattern=s)
    class(s) <- "summary.proximus"
    s
}

print.summary.proximus <- function(x, ...) {
    cat("approximates",x$nr,"x",x$nc, "matrix\n")
    cat("total Error:",format(x$error, digits=2), "\n")
    cat("total Fnorm:",format(x$fnorm, digits=2), "\n")
    cat("total  Jsim:",format(x$jsim,  digits=2), "\n")
    cat("total Valid:",x$valid,"\n")
    cat("Pattern Summary:\n")
    print(x$pattern[order(x$pattern$Size, decreasing=TRUE),], digits=2)
    invisible(x)
}

###

# Generate a matrix containing blocks of (overlapping) uniform 
# binary patterns on a noisy background. The perfect switch allows 
# for overlap between the first and last pattern block, making the 
# test case balanced.
#
# ceeboo 2005

rlbmat <- function(npat=4, rows=20, cols=12, over=4, noise=0.01, 
		   prob=0.8, perfect=FALSE) {
    
    rlmat <- function(nrow, ncol, prob=0.5) {
        x <- matrix(as.logical(runif(nrow*ncol) <= prob), ncol=ncol)
	    x
    }
    
    nrow <- npat * rows
    ncol <- cols * npat + over

    x <- rlmat(nrow, ncol, noise)

    r <- c <- 1
    while (r < nrow) {
        x[r:(r+rows-1), c:(c+cols+over-1)] <-
	        rlmat(rows, cols+over, prob)
	    r <- r + rows
	    c <- c + cols
    }
    # overlap first and last block, too
    if (perfect)
       x[(r-rows):(r-1), 1:over] <-
            rlmat(rows, cols, prob)
    x
}

###
