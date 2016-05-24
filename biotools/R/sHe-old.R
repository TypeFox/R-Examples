# S3method for "sHe"
sHe <- 
function(x, coord.cols = 1:2, marker.cols = 3:4,
	marker.type = c("codominant", "dominant"), 
	grid = NULL, latlong2km = TRUE, radius, nmin = NULL) 
   UseMethod("sHe")

# ------------------------------------------------
# default (data.frame or matrix)
sHe.default <- 
function(x, coord.cols = 1:2, marker.cols = 3:4,
	marker.type = c("codominant", "dominant"), 
	grid = NULL, latlong2km = TRUE, radius, nmin = NULL) 
{
    if (!inherits(x, c("data.frame", "matrix")))
       stop("'x' must be a data.frame o matrix")
    stopifnot(is.integer(coord.cols))
    stopifnot(is.integer(marker.cols))
    marker.type <- match.arg(marker.type)
    nomes <- rownames(x) <- as.character(rownames(x))

    # check coordinate system
    loc. <- as.matrix(x[, coord.cols])
    if (latlong2km) {
       loc <- as.matrix(latlong2grid(loc.))
    } else {
       loc <- loc.
    }
    d <- as.matrix(dist(loc))

    # check sampling points inside the prediction grid... or build a grid
    if (!is.null(grid)) {
       grid. <- as.matrix(grid)
       minmax <- apply(grid., 2, range)
       if ( any(loc.[, 1] < minmax[1, 1]) || 
          any(loc.[, 1] > minmax[2, 1]) ||
          any(loc.[, 2] < minmax[1, 2]) ||
          any(loc.[, 2] > minmax[2, 2]) ) {
             warning("some sampling coordinates seem to be outside the prediction grid")
       }
    } else {
       minmax <- apply(loc., 2, range)
       grid. <- as.matrix( expand.grid(x = seq(minmax[1, 1], minmax[2, 1], length.out = 50), 
          y = seq(minmax[1, 2], minmax[2, 2], length.out = 50)) )
    }
    if (latlong2km) {
       grid <- as.matrix(latlong2grid(grid.))
    } else {
       grid <- grid.
    }

    # selecting points within each centred coord of a grid
    mdis <- matrix(nrow = nrow(loc), ncol = nrow(grid))
    id <- list()
    for(i in 1:ncol(mdis)) {
       mdis[, i] <- sqrt(apply((loc - 
	    matrix(grid[i, ], nrow = nrow(loc), ncol = ncol(loc), 
             byrow = TRUE))^2, 1, sum))
       id[[i]] <- which(mdis[, i] <= radius)
    }
    n <- sapply(id, length)

    # organizing marker data
    markers <- as.matrix(x[, marker.cols])

    # progress bar... for the loop structure
    pb <- tkProgressBar(title = "Spatial Gene Diversity", 
       label = "CALCULATION PROGRESS", min = 0, max = length(id), width = 400L)

    # if codominant markers ------------------------------------------
    if (marker.type == "codominant") {
       if (ncol(markers) %% 2 != 0)
          stop("number of marker columns is not even!")
       o <- seq(1, ncol(markers), by = 2)

       # elements for loop and output
       MaxDist <- NULL
       mHe <- size <- bias <- matrix(nrow = length(id), ncol = length(o))

       # loop and calculation of gene diversity (He) in each grid point
       for (i in 1:length(id)) {
          MaxDist[i] <- ifelse(n[i] > 1, max(d[ id[[i]], id[[i]] ]), 0)
          for(j in 1:length(o)) {
          if (n[i] == 0) { 
             size[i, j] <- bias[i, j] <- mHe[i, j] <- 0
          } else if (!is.null(nmin) && n[i] < nmin) {
             size[i, j] <- bias[i, j] <- mHe[i, j] <- 0
          } else {
             size[i, j] <- sum(as.vector(markers[id[[i]], o[j]:(o[j]+1)]) > 0) / 2
             bias[i, j] <- 2*size[i, j] / ( 2*size[i, j] - 1 )
             mHe[i, j] <- ( 1 - sum( (table(as.vector(markers[id[[i]], o[j]:(o[j]+1)])) / 
                (2*size[i, j]) )^2) ) * bias[i, j]
          }
       }
       setTkProgressBar(pb, i, label = sprintf("CALCULATION PROGRESS (%.0f%%)", 
                100 * i/length(id)))
    }
    mHe[mHe < 0] <- mHe[is.na(mHe)] <- 0
    } else {
    # if dominant markers ---------------------------------------------
       # elements for loop and output
       MaxDist <- NULL
       mHe <- matrix(nrow = length(id), ncol = ncol(markers))

       # loop and calculation of gene diversity (He) in each grid point
       for (i in 1:length(id)) {
          MaxDist[i] <- ifelse(n[i] > 1, max(d[ id[[i]], id[[i]] ]), 0)
          for (j in 1:ncol(markers)) {
             if (n[i] == 0) { 
                mHe[i, j] <- 0
             } else if (!is.null(nmin) && n[i] < nmin) {
                mHe[i, j] <- 0
             } else {
                mHe[i, j] <- fHe(markers[id[[i]], j])
             }
          }       
          setTkProgressBar(pb, i, 
             label = sprintf("CALCULATION PROGRESS (%.0f%%)", 
             100 * i/length(id))) 
       }
    }

    # output
    uHe <- apply(mHe, 1, mean)
    SE <- apply(mHe, 1, sd) / sqrt(ncol(mHe))
    diversity <- data.frame(n, MaxDist, uHe, SE)
    out <- list(diversity = diversity, 
       mHe = mHe, grid = grid.)
    class(out) <- "sHe"
    Sys.sleep(0.5)
    close(pb)
    return(out)
}

# ------------------------------------------------
# aux function... uHe for COdominant markers
fHeco <- function(x)  # x is a two-columns matrix
{
   x <- as.matrix(x)
   n <- sum(!is.na(x)) / 2
   bias <- 2*n / (2*n - 1)
   He <- 1 - sum( (table(x) / (2*n))^2 )
   uHe <- He * bias
   return(uHe)
}

# ------------------------------------------------
# aux function... uHe for dominant markers
fHe <- function(x) { 
   x <- x[complete.cases(x)]
   n <- length(x)
   bias <- 2*n / (2*n - 1)
   fail <- sqrt(1 - mean(x))
   succ <- 1 - fail
   uHe <- 2 * succ * fail * bias
   return(uHe)
}

# ------------------------------------------------
# print method
print.sHe <- function(x, ...)
{
   print(summary(x$diversity), ...)
   invisible(x)
}

# ------------------------------------------------
# plot method (lattice::levelplot)
plot.sHe <- function(x, ...)
{
   levelplot(x$diversity$uHe ~ x$grid[, 1] * x$grid[, 2],
      col.regions = rev(heat.colors(100)),
      main = "Gene Diversity Heat Map", ...)
}
