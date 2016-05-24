###########################################################################
##                                                                       ##
## distance - function to compute distances between samples              ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 13-Oct-2007                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
##                                                                       ##
###########################################################################
## x = training data, y = fossil data
oldDistance <- function(x, ...) UseMethod("oldDistance")

oldDistance.join <- function(x, ...)
  {
    if(!inherits(x, "join"))
      stop("This method should only be used on objects of class 'join'")
    if(inherits(x, "data.frame")) {
      distance.default(x, ...)
    } else {
      if(length(x) != 2)
        warning("Object contains more than 2 data sets.\n  Only the first 2 data sets used")
      distance.default(x[[1]], x[[2]], ...)
    }
  }

oldDistance.default <- function(x, y,
                             method = c("euclidean", "SQeuclidean", "chord",
                               "SQchord", "bray", "chi.square", "SQchi.square",
                               "information", "chi.distance", "manhattan",
                               "kendall", "gower", "alt.gower", "mixed"),
                             fast = TRUE,
                             weights = NULL, R = NULL, ...)
  {
    euclidean <- function(x, y)
      {
        sqrt(sum((x - y)^2))
      }
    SQeuclidean <- function(x, y)
      {
        sum((x - y)^2)
      }
    chord <- function(x, y)
      {
        x <- sqrt(x); y <- sqrt(y)
        euclidean(x, y)
      }
    SQchord <- function(x, y)
      {
        x <- sqrt(x); y <- sqrt(y)
        SQeuclidean(x, y)
      }
    bray <- function(x, y)
      {
        sum(abs(x - y)) / sum(x + y)
      }
    chi.square <- function(x, y)
      {
        inds <- !(x == 0 & y == 0)
        sqrt(sum(((x[inds] - y[inds])^2) / (x[inds] + y[inds])))
      }
    SQchi.square <- function(x, y)
      {
        inds <- !(x == 0 & y == 0)
        sum(((x[inds] - y[inds])^2) / (x[inds] + y[inds]))
      }
    information <- function(x, y)
      {
        XY <- x + y
        A <- x * log2((2 * x) / XY)
        B <- y * log2((2 * y) / XY)
        sum(A, na.rm = TRUE) + sum(B, na.rm = TRUE)
      }
    chi.distance <- function(x, y, colsum)
      {
        sqrt(sum(((x - y)^2) / (colsum / sum(colsum))))
      }
    manhattan <- function(x, y)
      {
        sum(abs(x - y))
      }
    kendall <- function(x, y, maxi)
      {
          ## the sum in the else isn't right if x == y
          ## then the dissimilarity should be 0
          out <- if (isTRUE(all.equal(sum(x-y), 0))) {
              0
          } else {
              sum(maxi - pmin(x, y))
          }
          out
      }
    gower <- function(x, y, maxi, mini)
      {
        sum(abs(x - y) / (maxi - mini), na.rm = TRUE)
      }
    alt.gower <- function(x, y, maxi, mini)
      {
        sqrt(2 * sum(abs(x - y) / (maxi - mini), na.rm = TRUE))
      }
    mixed <- function(x, y, facs, weights, R)
      {
        weights[is.na(x)] <- 0
        weights[is.na(y)] <- 0
        if(any(!facs)) {
          quant <- 1 - (abs(x[!facs] - y[!facs]) / R[!facs])
          quant <- quant * weights[!facs]
        }
        if(any(facs)) {
          factors <- ifelse(x[facs] == y[facs], 1, 0)
          factors <- factors * weights[facs]
        }
        if(any(!facs)) {
          if(any(facs))
            retval <- sum(factors, quant, na.rm = TRUE) / sum(weights)
          else
            retval <- sum(quant, na.rm = TRUE) / sum(weights)
        } else {
          retval <- sum(factors, na.rm = TRUE) / sum(weights)
        }
        return(1- retval)
      }
    Dist <- function(y, x, method, ...)#colsum = NULL)
      {
        dotargs <- list(...)
        switch(method,
               euclidean = apply(x, 1, euclidean, y),
               SQeuclidean = apply(x, 1, SQeuclidean, y),
               chord = apply(x, 1, chord, y),
               SQchord = apply(x, 1, SQchord, y),
               bray = apply(x, 1, bray, y),
               chi.square = apply(x, 1, chi.square, y),
               SQchi.square = apply(x, 1, SQchi.square, y),
               information = apply(x, 1, information, y),
               chi.distance = apply(x, 1, chi.distance, y,
                 colsum = dotargs$colsum),
               manhattan = apply(x, 1, manhattan, y),
               kendall = apply(x, 1, kendall, y, maxi = dotargs$maxi),
               gower = apply(x, 1, gower, y, maxi = dotargs$maxi,
                 mini = dotargs$mini),
               alt.gower = apply(x, 1, alt.gower, y, maxi = dotargs$maxi,
                 mini = dotargs$mini),
               mixed = apply(x, 1, mixed, y, facs = dotargs$facs,
                 weights = dotargs$weights, R = dotargs$R)
               )
      }
    if(missing(method))
        method <- "euclidean"
    method <- match.arg(method)
    y.miss <- FALSE
    if(missing(y)) {
        y.miss <- TRUE
        y <- x
    }
    n.vars <- ncol(x)
    if(method == "mixed") {
        ## are same columns in x and y factors
        facs.x <- sapply(as.data.frame(x), is.factor, USE.NAMES = FALSE)
        facs.y <- sapply(as.data.frame(y), is.factor, USE.NAMES = FALSE)
        if(sum(facs.x - facs.y) > 0) {
            stop("Different columns (species) are coded as factors in 'x' and 'y'")
            ## levels of factors also need to be the same
            for(i in seq_along(facs.x)[facs.x]){
                if(!identical(levels(x[,i]), levels(y[,i])))
                    stop("The levels of one or more factors in 'x' and 'y'\ndo not match.\nConsider using 'join(x, y)'. See '?join'")
            }
        }
    } else {
        ## we do this even if no y as it is harmless
        facs.x <- facs.y <- rep(FALSE, n.vars)
    }
    x.names <- rownames(x)
    x <- data.matrix(x)
    y.names <- rownames(y)
    y <- data.matrix(y)
    ## Do we want to remove NAs? Yes if gower, alt.gower and mixed,
    ## but fail for others
    NA.RM <- FALSE
    if(method %in% c("gower", "alt.gower", "mixed"))
        NA.RM <- TRUE
    ## check if any empty species, drop them
    colsumx <- colSums(x, na.rm = NA.RM)
    colsumy <- colSums(y, na.rm = NA.RM)
    ## NO - this causes problems if you merge data
    if (any(DROP <- (colsumx <= 0 & colsumy <= 0) & !facs.x)) {
        ##x <- x[, (colsumx > 0 | colsumy > 0) | facs.x, drop = FALSE]
        ##y <- y[, (colsumx > 0 | colsumy > 0) | facs.x, drop = FALSE]
        ##warning("Some species contain no data and were removed from data matrices.\n")
    }
    if(method == "chi.distance")
        colsum <- colSums(join(as.data.frame(x),as.data.frame(y), split = FALSE))
    if(method == "mixed") {
      ## sort out the weights used, eg the Kroneker's Deltas
      ## weights must be NULL or numeric vector of length == ncol(x)
      if(is.null(weights))
        weights <- rep(1, n.vars)
      else {
        if(length(weights) != n.vars)
          stop("'weights' must be of length 'ncol(x)'")
      }
    }
    dimx <- dim(x)
    dimy <- dim(y)
    if(method %in% c("gower", "alt.gower", "mixed", "kendall")) {
        ## for these methods, need max for each var and mins if not
        ## doing Kendall
        NA.RM <- TRUE
        if(method == "kendall")
            NA.RM <- FALSE
        ## need to account for a single site (matrix with 1 row)
        ## that might have NA. Specifically not allowed in Kendall
        ## but OK in the other methods handled here
        maxX <- if(any(apply(is.na(x), 2, all))) {
            x
        } else {
            apply(x, 2, max, na.rm = NA.RM)
        }
        maxY <- if(any(apply(is.na(y), 2, all))) {
            y
        } else {
            apply(y, 2, max, na.rm = NA.RM)
        }
        ##maxi <- apply(rbind(maxX, maxY), 2, max, na.rm = NA.RM)
        maxi <- pmax(maxX, maxY)
        if(method %in% c("gower", "alt.gower", "mixed")) {
            ## need the mins of each variable
            ## need to account for a single site (matrix with 1 row)
            ## that might have NA. Specifically not allowed in Kendall
            ## but OK in the other methods handled here
            minX <- if(any(apply(is.na(x), 2, all))) {
                x
            } else {
                apply(x, 2, min, na.rm = NA.RM)
            }
            minY <- if(any(apply(is.na(y), 2, all))) {
                y
            } else {
                apply(y, 2, min, na.rm = NA.RM)
            }
            mini <- apply(rbind(minX, minY), 2, min, na.rm = NA.RM)
            ## compute R - the range - if not supplied
            ## if R supplied then validate
            if(is.null(R)) {
                R <- maxi - mini
            } else {
                if(length(R) != n.vars)
                    stop("'R' must be of length 'ncol(x)'")
            }
        }
    }
    dimnames(x) <- dimnames(y) <- NULL
    if(method == "chi.distance") {
      y <- y / rowSums(y)
      x <- x / rowSums(x)
      res <- apply(y, 1, Dist, x, method, colsum = colsum)
    } else if (method == "kendall") {
      res <- apply(y, 1, Dist, x, method, maxi = maxi)
    } else if (method %in% c("gower", "alt.gower")) {
      res <- apply(y, 1, Dist, x, method, maxi = maxi, mini = mini)
    } else if (method == "mixed") {
      res <- apply(y, 1, Dist, x, method, facs = facs.x, weights = weights,
                   R = R)
    } else if (method %in% c("euclidean","SQeuclidean","chord","SQchord","bray") &&
               fast == TRUE && y.miss == TRUE) {
      if(method == "bray") {
        res <- as.matrix(vegdist(x, method = "bray"))
      } else {
        if(method %in% c("chord", "SQchord"))
          res <- as.matrix(dist(sqrt(x), method = "euclidean"))
        else
          res <- as.matrix(dist(x, method = "euclidean"))
        if(method %in% c("SQeuclidean", "SQchord"))
          res <- res^2
      }
    } else {
      res <- apply(y, 1, Dist, x, method)
    }
    if(is.null(dim(res))) {
      names(res) <- x.names
    } else {
      colnames(res) <- y.names
      rownames(res) <- x.names
    }
    attr(res, "method") <- method
    attr(res, "type") <- if(y.miss) "symmetric" else "asymmetric"
    class(res) <- c("distance","matrix")
    return(res)
  }
