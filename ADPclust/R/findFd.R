#' @title Find f(x) and delta(x) from a Data Set
#' 
#' @description This function finds f(x) and delta(x) of each point x, where f(x) is 
#' the density estimation at each point, and delta(x) measures the 
#' distance between x to the closest point y such that f(y) > f(x). 
#' The output is used in 'findCluster()' and 'findClusterAuto()'.
#' 
#' @param dat a numeric matrix where rows are samples and columns are features.
#' @param h nonnegative number or NULL (default). Bandwidth in density estimation. If NULL h is automatically selected.
#' @param htype integer vector. Number of clusters. Default = 5.
#' @param dmethod character string describing distance measures to use in dist(.. method = dmethod). 
#'
#' @param fdelta character string that specifies the method to estimate densities at each data point. Default is "mnorm" stands for multivariate Gaussian density estimation. Other options include
#' @param verbose if TRUE progress will be displayed.#'   
#' @return an R object with class "rd". If length(h) = 1 then "rd" is a list of 
#' the following items:
#' \itemize{
#' \item f: vector of f's
#' \item delta: vector of delta's
#' \item dat: the original data matrix
#' \item distm: distance matrix calculated from data (Euclidean distance))
#' \item dc: bandwith dc.
#' }
#' If h == NULL, then the function returns a list of length 15, where each 
#' item is a list corresponding to one h value. 
#' 
#' @examples
#' ## data(clust3)
#' ##findFd(clust3, h = 19, plot = TRUE)

findFd <- function(dat, h = NULL,
                   dmethod = "euclidean",
                   htype = "AMISE",
                   fdelta = "mnorm",
                   verbose = FALSE)
{
    n <- nrow(dat)
    p <- ncol(dat)

    ## Calcuate the distance matrix from standarized data.
    if(verbose) cat("Calculating distm. ")
    if(fdelta == "mnorm"){
        ## distm.std <- as.matrix(dist(scale(dat, center = FALSE, scale = TRUE),
        ##                             method = "euclidean"))
        sds <- apply(dat, 2, sd)
        distm <- as.matrix(dist(scale(dat, center = FALSE, scale = sds),
                                    method = "euclidean", upper = TRUE))
    }else{
        distm <- as.matrix(dist(dat, upper = TRUE, method = dmethod))
    }
    if(verbose) cat("Done.\n")
    
    find_fd <- function(h){
        ## Find f
        if(verbose) cat("Calculating f for h = ", h, "... ")
        if(fdelta == "normkernel")
            f <- 1/(h * sqrt(2 * pi)) * rowSums(exp(-(distm/h)^2/2))
        else if(fdelta == "weighted")
            f <- rowSums(exp(-(distm/h)^2))
        else if(fdelta == "count")
            f <- rowSums(distm < h) - 1
        else if(fdelta == "mnorm"){
            f <- rowSums(exp(-(distm / h) ^ 2 / 2))
            if(all(f == 1)) stop("All f = 1, the bandwidth h =", h, " is too small.")
        }else{
            stop("Wrong fdelta, try 'normkernel', 'weighted', 'count' or 'mnorm' (recommended).")
        }
        if(verbose) cat("Done.\n")
        ## Find delta
        if(verbose) cat("Calculating delta for h = ", h, "... ")
        if(fdelta == "count"){
            f1 <- rank(f, ties.method = "first") # Break ties in f
            delta <- apply(distm / outer(f1, f1, FUN = ">"), 2, min, na.rm = TRUE)
            loc.max <- which.max(delta)
            delta[loc.max] <- max(delta[-loc.max]) # Equation in the Matlab cqode
        }else if(fdelta == "mnorm"){
            time1 <- Sys.time()
            f.order <- order(f, decreasing = TRUE)
            delta <- rep(NA, n)
            delta[f.order[1]] <- Inf
            for(i in 2:length(f.order)){
                delta[f.order[i]] <- min(distm[f.order[i], f.order[1:(i - 1)]])
            }
            delta[f.order[1]] <- max(delta[-f.order[1]])
            
            ## time2 <- Sys.time()
            ## frank <- rank(f, ties.method = "random")
            ## delta <- rep(NA, n)
            ## for(i in 1:n){
            ##     delta[i] <- min(distm.std[i, frank > frank[i]])
            ## }
            
            ## time3 <- Sys.time()
            ## delta <- apply(distm.std / outer(f, f, FUN = ">"), 2, min, na.rm = TRUE)
            
            ## time4 <- Sys.time()
            ## cat("loops uses", time2 - time1,
            ##     "loops2 uses", time3 - time2,
            ##     "apply uses", time4 - time3, "\n")
        }else{
            delta <- apply(distm / outer(f, f, FUN = ">"), 2, min, na.rm = TRUE)
            loc.max <- which.max(delta)
            delta[loc.max] <- max(delta[-loc.max]) # Equation in the Matlab cqode
        }
        if(verbose) cat("Done.\n")
        return(list(f = f, delta = delta))
    }

    ##--------------------------
    ## Main function
    ##--------------------------
    if(is.null(h) | length(h) > 1){
        if(is.null(h)){
            if(fdelta != "mnorm")
                stop("h must be a nonnegative number.")
            else{
                browser()
                if(htype == "AMISE")
                    h <- AMISE(p, n)
                else if(htype == "ROT")
                    h <- ROT(p, n)
                else
                    stop("h is missing. Wrong htype. Give h, or set htype = 'AMISE' or htype = 'ROT'.")
            }
            h <- seq(h / 3, 3 * h, length.out = 10)
        }else{
            if(!is.numeric(h) || !all(h > 0)) stop("h must be nonnegative number(s).")
        }
        fds <- lapply(h, find_fd)
        res <- list(fd = fds, dat = dat, distm = distm, h = h)
    }else{
        if(!is.numeric(h) || h <= 0) stop("h must be nonnegative number(s)")
        ## h known
        fd <- list(find_fd(h))
        res <- list(fd = fd, dat = dat, distm = distm, h = h)
    }
    class(res) <- c("fd", "list")
    attributes(res)$state <- ifelse(is.null(h) | length(h) > 1, "single.h", "multiple.h")
    invisible(res)
}
