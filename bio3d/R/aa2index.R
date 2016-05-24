"aa2index" <-
function (aa, index = "KYTJ820101", window = 1) {
    if (!is.vector(aa)) 
        stop("aa2index: non vector argument")
    if (!is.numeric(window) || window < 1) 
        stop("aa2index: 'window' must be numeric and positive")
    if (window >= length(aa)) 
        stop("aa2index: 'window' must be smaller than the sequence length")
     # Use LazyData to import data - changed Jul 23, 2013
#    if (!exists("aa.index")) 
#        data(aa.index)
    aa.index = bio3d::aa.index
    if (is.numeric(index)) {
        if (index > length(names(aa.index))) {
            stop("aa2index: 'index' number does not exist")
        }
    } else {
        if (!is.element(index, names(aa.index))) {
            stop("aa2index: 'index' name does not exist")
        }
    }
    x <- aa.index[[index]]$I[aa]
    if (window == 1) {
        y <- x
    } else {
        n <- length(x)
        y <- rep(NA, n)
        w <- ceiling(window/2)
        
        if ( (window %% 2) == 0 ) {
          from <- w
          to <- n - w
          y[from:to] <- sapply(from:to, function(i)
                               mean(x[(i - (w-1)):(i+w)], na.rm=TRUE))
          if (from-1 > 0) {
            y[1:(from-1)] <- sapply(1:(from-1), function(i)
                                    mean(x[1: (i + w)], na.rm=TRUE))
          }
          y[(to+1):n] <- sapply((to+1):n, function(i)
                                mean(x[(i- (w-1)): n], na.rm=TRUE))
        } else {
          from <- w
          to <- n - (w-1)
          y[from:to] <- sapply(from:to, function(i)
                               mean(x[(i - (w-1)):(i + (w-1))], na.rm=TRUE))
          y[1:(from-1)] <- sapply(1:(from-1), function(i)
                                  mean(x[1: (i + (w-1))], na.rm=TRUE))   
          y[(to+1):n] <- sapply((to+1):n, function(i)
                                mean(x[(i - (w-1)): n], na.rm=TRUE))
        }
        y <- round(y,2)
        names(y) <- aa
    }
    return(y)
}


