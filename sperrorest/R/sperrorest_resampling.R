
#' Partition the data for a (non-spatial) cross-validation
#'
#' \code{partition.cv} creates a \code{\link{represampling}} object for \code{length(repetition)}-repeated \code{nfold}-fold cross-validation.
#' @name partition.cv
#' @aliases partition.cv
#' @param data \code{data.frame} containing at least the columns specified by \code{coords}
#' @param coords (ignored by \code{partition.cv})
#' @param nfold number of partitions (folds) in \code{nfold}-fold cross-validation partitioning
#' @param repetition numeric vector: cross-validation repetitions to be generated. Note that this is 
#'   not the number of repetitions, but the indices of these repetitions. E.g., use 
#'   \code{repetition=c(1:100)} to obtain (the 'first') 100 repetitions, and 
#'   \code{repetition=c(101:200)} to obtain a different set of 100 repetitions.
#' @param seed1 \code{seed1+i} is the random seed that will be used by \code{\link{set.seed}} 
#'   in repetition \code{i} (\code{i} in \code{repetition}) to initialize the random number 
#'   generator before sampling from the data set.
#' @param return.factor if \code{FALSE} (default), return a \code{\link{represampling}} object; if \code{TRUE} (used internally by other \code{sperrorest} functions), return a \code{list} containing factor vectors (see Value)
#' @details This function does not actually perform a cross-validation or partition the data set itself; 
#' it simply creates a data structure containing the indices of training and test samples.
#' @return If \code{return.factor=FALSE} (the default), a \code{\link{represampling}} object. 
#' Specifically, this is a (named) list of \code{length(repetition)} \code{resampling} objects.
#' Each of these \code{\link{resampling}} objects is a list of length \code{nfold} corresponding to the folds.  
#' Each fold is represented by a list of containing the components \code{train} and \code{test},
#' specifying the indices of training and test samples (row indices for \code{data}).
#' If \code{return.factor=TRUE} (mainly used internally), a (named) list of length \code{length(repetition)}. 
#' Each component of this list is a vector of length \code{nrow(data)} of type \code{factor}, specifying 
#' for each sample the fold to which it belongs. The factor levels are \code{factor(1:nfold)}.
#' @seealso \code{\link{sperrorest}}, \code{\link{represampling}}
#' @examples
#' data(ecuador)
#' ## non-spatial cross-validation:
#' resamp = partition.cv(ecuador, nfold = 5, repetition = 1:2)
#' plot(resamp, ecuador)
#' # first repetition, second fold, test set indices:
#' idx = resamp[["1"]][[2]]$test
#' # test sample used in this particular repetition and fold:
#' ecuador[ idx , ]
#' @export
partition.cv = function(data, coords = c("x","y"), nfold = 10, repetition = 1, seed1 = NULL,
    return.factor = FALSE)
{  
    resampling = list()

    for (cnt in repetition) {
        if (!is.null(seed1)) set.seed(seed1 + cnt)
        resampler = sample(rep(sample(nfold), length = nrow(data)), size = nrow(data))
        resampler = factor(resampler)
        if (!return.factor) resampler = as.resampling(resampler)
        resampling[[as.character(cnt)]] = resampler
    }
    if (!return.factor) resampling = as.represampling(resampling)

    return(resampling)
}


#' Partition the data for a stratified (non-spatial) cross-validation
#'
#' \code{partition.cv.strat} creates a set of sample indices corresponding to cross-validation test and training sets.
#' @name partition.cv.strat
#' @inheritParams partition.cv
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param strat character: column in \code{data} containing a factor variable over which the partitioning should be stratified; or factor vector of length \code{nrow(data)}: variable over which to stratify
#' @return A \code{\link{represampling}} object, see also \code{\link{partition.cv}}. \code{partition.strat.cv}, 
#' however, stratified with respect to the variable \code{data[,strat]}; i.e., cross-validation partitioning 
#' is done within each set \code{data[data[,strat]==i,]} (\code{i} in \code{levels(data[,strat])}), and the 
#' \code{i}th folds of all levels are combined into one cross-validation fold.
#' @seealso \code{\link{sperrorest}}, \code{\link{as.resampling}}, \code{\link{resample.strat.uniform}}
#' @examples
#' data(ecuador)
#' parti = partition.cv.strat(ecuador, strat = "slides", nfold = 5, repetition = 1)
#' idx = parti[["1"]][[1]]$train
#' mean(ecuador$slides[idx]=="TRUE") / mean(ecuador$slides=="TRUE")
#' # always == 1
#' # Non-stratified cross-validation:
#' parti = partition.cv(ecuador, nfold = 5, repetition = 1)
#' idx = parti[["1"]][[1]]$train
#' mean(ecuador$slides[idx]=="TRUE") / mean(ecuador$slides=="TRUE")
#' # close to 1 because of large sample size, but with some random variation
#' @export
partition.cv.strat = function(data, coords = c("x","y"), nfold = 10,
    return.factor = FALSE, repetition = 1, seed1 = NULL, strat)
{  
    repres = list()

    stopifnot( (length(strat) == 1) | (length(strat) == nrow(data)) )
    if (length(strat) == 1) strat = data[,strat]
    stopifnot(is.factor(strat))
    # Can't split into nfold partitions if there are less than nfold samples within a stratum:
    minstrat = min(tapply(strat, strat, length))
    stopifnot(minstrat >= nfold)

    for (cnt in repetition) {    
        if (!is.null(seed1)) set.seed(seed1 + cnt)
        fac = rep(NA, nrow(data))
        for (lev in levels(strat)) {
            nstrat = sum( sel <- (strat == lev))
            fac[sel] = sample(rep(sample(nfold), length = nstrat), size = nstrat)
        }
        fac = factor(fac)
        if (!return.factor) fac = as.resampling(fac)
        repres[[as.character(cnt)]] = fac
    }
    if (!return.factor) repres = as.represampling(repres)

    return(repres)
}

#' Partition the data for a (non-spatial) leave-one-factor-out cross-validation based on a given, fixed partitioning
#'
#' \code{partition.factor} creates a \code{\link{represampling}} object, i.e. a set of sample indices defining cross-validation test and training sets.
#' @inheritParams partition.cv
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param fac either the name of a variable (column) in \code{data}, or a vector of type factor and length \code{nrow(data)} that contains the partitions to be used for defining training and test samples
#' @return A \code{\link{represampling}} object, see also \code{\link{partition.cv}} for details.
#' @note In this partitioning approach, all \code{repetition}s are identical and therefore pseudo-replications.
#' @seealso \code{\link{sperrorest}}, \code{\link{partition.cv}}, \code{\link{as.resampling.factor}}
#' @examples 
#' data(ecuador)
#' # I don't recommend using this partitioning for cross-validation,
#' # this is only for demonstration purposes:
#' breaks = quantile(ecuador$dem, seq(0,1,length=6))
#' ecuador$zclass = cut(ecuador$dem, breaks, include.lowest=TRUE)
#' summary(ecuador$zclass)
#' parti = partition.factor(ecuador, fac = "zclass")
#' plot(parti,ecuador)
#' summary(parti)
#' @export
partition.factor = function(data, coords = c("x", "y"), fac, return.factor = FALSE,
    repetition = 1)
{
    if (length(fac) == 1 && is.character(fac)) fac = data[,fac]
    fac = factor(fac)
    if (!return.factor) fac = as.resampling(fac)
    represmp = list()
    for (cnt in repetition)
        represmp[[as.character(cnt)]] = fac
    if (!return.factor) represmp = as.represampling(represmp)
    return(represmp)
}

#' Partition the study area into rectangular tiles
#'
#' \code{partition.tiles} divides the study area into a specified number of rectangular tiles. Optionally small partitions can be
#' merged with adjacent tiles to achieve a minimum number or percentage of samples in each tile.
#' @inheritParams partition.cv
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param dsplit optional vector of length 2: equidistance of splits in (possibly rotated) x direction (\code{dsplit[1]}) and y direction (\code{dsplit[2]}) used to define tiles. If \code{dsplit} is of length 1, its value is recycled. Either \code{dsplit} or \code{nsplit} must be specified.
#' @param nsplit optional vector of length 2: number of splits in (possibly rotated) x direction (\code{nsplit[1]}) and y direction (\code{nsplit[2]}) used to define tiles. If \code{nsplit} is of length 1, its value is recycled.
#' @param rotation indicates whether and how the rectangular grid should be rotated; random rotation is only between \code{-45} and \code{+45} degrees
#' @param user.rotation if \code{rotation="user"}, angles (in degrees) by which the rectangular grid is to be rotated in each repetition. Either a vector of same length as \code{repetition}, or a single number that will be replicated \code{length(repetition)} times.
#' @param offset indicates whether and how the rectangular grid should be shifted by an offset
#' @param user.offset if \code{offset="user"}, a list (or vector) of two components specifying a shift of the rectangular grid in (possibly rotated) x and y direction. The offset values are relative values, a value of \code{0.5} resulting in a one-half tile shift towards the left, or upward. If this is a list, its first (second) component refers to the rotated x (y) direction, and both components must have same length as \code{repetition} (or length 1). If a vector of length 2 (or list components have length 1), the two values will be interpreted as relative shifts in (rotated) x and y direction, respectively, and will therefore be recycled as needed (\code{length(repetition)} times each).
#' @param reassign logical (default \code{TRUE}): if \code{TRUE}, 'small' tiles (as per \code{min.frac} and \code{min.n} arguments and \code{\link{get.small.tiles}}) are merged with (smallest) adjacent tiles. If \code{FALSE}, small tiles are 'eliminated', i.e. set to \code{NA}.
#' @param min.frac numeric >=0, <1: minimum relative size of partition as percentage of sample; argument passed to \code{\link{get.small.tiles}}. Will be ignored if \code{NULL}.
#' @param min.n integer >=0: minimum number of samples per partition; argument passed to \code{\link{get.small.tiles}}. Will be ignored if \code{NULL}.
#' @param iterate argument to be passed to \code{\link{tile.neighbors}}
#' @return A \code{\link{represampling}} object. Contains \code{length(repetition)} \code{\link{resampling}} 
#' objects as repetitions. The exact number of folds / test-set tiles within each \code{\link{resampling}} 
#' objects depends on the spatial configuration of the data set and possible cleaning steps 
#" (see \code{min.frac}, \code{min.n}).
#' @note Default parameter settings may change in future releases. This function, especially the rotation and shifting part of it and the algorithm for cleaning up small tiles is still a bit experimental. Use with caution.
#' For non-zero offsets (\code{offset!="none")}), the number of tiles may actually be greater than \code{nsplit[1]*nsplit[2]} because of fractional tiles lurking into the study region. \code{reassign=TRUE} with suitable thresholds is therefore recommended for non-zero (including random) offsets.
#' @seealso \code{\link{sperrorest}}, \code{\link{as.resampling.factor}}, \code{\link{get.small.tiles}}, \code{\link{tile.neighbors}}
#' @examples
#' data(ecuador)
#' parti = partition.tiles(ecuador, nsplit = c(4,3), reassign = FALSE)
#' plot(parti,ecuador)
#' summary(parti) # tile A4 has only 55 samples
#' # same partitioning, but now merge tiles with less than 100 samples to adjacent tiles:
#' parti2 = partition.tiles(ecuador, nsplit = c(4,3), reassign = TRUE, min.n = 100)
#' plot(parti2,ecuador)
#' summary(parti2)
#' # tile B4 (in 'parti') was smaller than A3, therefore A4 was merged with B4, not with A3
#' # now with random rotation and offset, and tiles of 2000 m length:
#' parti3 = partition.tiles(ecuador, dsplit = 2000, offset = "random", rotation = "random", reassign = TRUE, min.n = 100)
#' plot(parti3,ecuador)
#' summary(parti3)
#' @export
partition.tiles = function(data, coords = c("x","y"), dsplit = NULL, nsplit = NULL,
    rotation = c("none","random","user"), user.rotation,
    offset = c("none","random","user"), user.offset,
    reassign = TRUE, min.frac = 0.025, min.n = 5, iterate = 1,
    return.factor = FALSE, repetition = 1, seed1 = NULL)
{
    # Some basic argument checks:
    stopifnot(is.numeric(min.frac) && length(min.frac)==1)
    stopifnot(is.numeric(min.n) && length(min.n)==1)
    stopifnot(is.numeric(iterate) && length(iterate)==1)
    stopifnot(!is.null(nsplit) | !is.null(dsplit))
    if (!is.null(nsplit)) {
        stopifnot(is.numeric(nsplit) && length(nsplit)<=2)
    } else stopifnot(is.numeric(dsplit) && length(dsplit)<=2)

    # Prepare rotation angles, if applicable:
    rotation = match.arg(rotation)
    stopifnot(xor(rotation == "user", missing(user.rotation)))
    if (rotation == "none") {
        phi = rep(0, length(repetition))
    } else if (rotation == "random") {
        phi = runif(-45, 45, n = length(repetition))
    } else if (rotation == "user") {
        if (length(user.rotation) == 1) user.rotation = rep(user.rotation, length(repetition))
        stopifnot(length(user.rotation) == length(repetition))
        phi = user.rotation
    }
    names(phi) = as.character(repetition)
    
    # This will make matrix multiplication (rotation) numerically better conditioned:
    data[,coords[1]] = data[,coords[1]] - mean(data[,coords[1]])
    data[,coords[2]] = data[,coords[2]] - mean(data[,coords[2]])

    offset = match.arg(offset)
    stopifnot(xor(offset == "user", missing(user.offset)))
    if (offset == "none") {
        x.shift = y.shift = rep(0, length(repetition))
    } else if (offset == "random") {
        x.shift = runif(0, 1, n = length(repetition))
        y.shift = runif(0, 1, n = length(repetition))
    } else if (offset == "user") {
        if (is.vector(user.offset) && length(user.offset) == 2)
            user.offset = list(user.offset[1], user.offset[2])
        stopifnot(is.list(user.offset) && length(user.offset) == 2)
        # Recycle offsets as needed:
        if (length(user.offset[[1]]) == 1) user.offset[[1]] = rep(user.offset[[1]], length(repetition))
        if (length(user.offset[[2]]) == 1) user.offset[[2]] = rep(user.offset[[2]], length(repetition))
        # Got enough user.offsets?
        stopifnot(length(user.offset[[1]]) == length(repetition))
        stopifnot(length(user.offset[[2]]) == length(repetition))
        # Valid range, [0,1]?
        stopifnot( min(user.offset[[1]] >= 0 & max(user.offset[[1]]) <= 1) )
        stopifnot( min(user.offset[[2]] >= 0 & max(user.offset[[2]]) <= 1) )
        x.shift = user.offset[[1]]
        y.shift = user.offset[[2]]
    }
    names(x.shift) = as.character(repetition)
    names(y.shift) = as.character(repetition)

    if (!is.null(nsplit))
        if (length(nsplit) == 1)  nsplit = c(nsplit, nsplit)
    if (!is.null(dsplit))
        if (length(dsplit) == 1)  dsplit = c(dsplit, dsplit)

    resampling = list()
    for (cnt in repetition) {

        if (!is.null(seed1)) set.seed(seed1 + cnt)

        # Prepare the arguments and data:

        if (rotation != "none") {        
            R = phi[as.character(cnt)] * 180 / pi
            R = matrix(c(cos(R), -sin(R), sin(R), cos(R)), ncol = 2)
            xy = R %*% t(data[,coords])
            x = xy[1,]
            y = xy[2,]
        } else {
            x = data[,coords[1]]
            y = data[,coords[2]]
        }
        x.range = range(x)
        y.range = range(y)
        
        if (!is.null(nsplit)) {
            x.delta = diff(x.range)/nsplit[1]
            y.delta = diff(y.range)/nsplit[2]
            my.nsplit = nsplit
        } else { # if !is.null(dsplit)
            x.delta = dsplit[1]
            y.delta = dsplit[2]
        }
        # Apply offsets:
        if (offset != "none") {
            # Widen the range and increase nsplit to allow for "lurking" tiles:
            x.range[2] = x.range[2] + x.delta
            y.range[2] = y.range[2] + y.delta
            x.range = x.range - x.delta * (x.shift[as.character(cnt)])
            y.range = y.range - y.delta * (y.shift[as.character(cnt)])
            if (is.null(dsplit)) my.nsplit = my.nsplit + 1
        }
        
        # Calculate x and y splits:
        if (is.null(dsplit)) {
            x.split = seq(x.range[1], x.range[2], length = my.nsplit[1] + 1)
            y.split = seq(y.range[1], y.range[2], length = my.nsplit[2] + 1)
        } else {
            x.split = seq(x.range[1], x.range[2] + x.delta, by = x.delta)
            y.split = seq(y.range[1], y.range[2] + y.delta, by = y.delta)
            my.nsplit = c(length(x.split)-1, length(y.split)-1)
        }
            
        # Group data into tiles, i.e. assign tile labels to samples:
        tile = rep(NA, nrow(data))
        for (ix in 1:my.nsplit[1]) {
            # Intervals are normally open to the left, except the first one:
            if (ix == 1) {
                sel.x = (x >= x.split[ix]) & (x <= x.split[ix+1])
            } else {
                sel.x = (x > x.split[ix]) & (x <= x.split[ix+1])
            }
            for (iy in 1:my.nsplit[2]) {
                if (iy == 1) {
                    sel.y = (y >= y.split[iy]) & (y <= y.split[iy+1])
                } else {
                    sel.y = (y > y.split[iy]) & (y <= y.split[iy+1])
                }
                # Assign tile name to samples:
                if (any( sel.x & sel.y ))
                    tile[ sel.x & sel.y ] = as.character(as.tilename(c(ix,iy)))
            }
        }
        tile = factor(tile)
    
        # Identify and process small tiles:
        s.tiles = get.small.tiles(tile, min.n = min.n, min.frac = min.frac)
        if (length(s.tiles) > 0) { # any small tiles?
            if (reassign) {
                # Merge small tiles with neighbors:
                ignore = c()
                # Repeat until no small tiles are left:
                while ((length(s.tiles) > 0) & (length(levels(tile)) > 1)) {
                    # Start with smallest small tile:
                    nbrs = tile.neighbors(s.tiles[1], tileset = levels(tile), iterate=iterate)
                    if (length(nbrs) == 0) {
                        ignore = c(ignore, as.character(s.tiles[1]))
                    } else {
                        # Merge tile with smallest neighbour to keep tile sizes balanced:
                        n.tile = tapply(tile, tile, length)
                        s.nbr = nbrs[ which.min(n.tile[nbrs]) ]
                        tile[ tile == s.tiles[1] ] = s.nbr
                        tile = factor(as.character(tile))
                    }
                    # Update small tiles list:
                    s.tiles = get.small.tiles(tile, min.n = min.n, min.frac = min.frac, ignore = ignore)
                }
            } else {
                # Just eliminate small tiles:
                tile[ tile %in% s.tiles ] = NA
                tile = factor(as.character(tile))
            }
        }
        
        if (!return.factor) tile = as.resampling(tile)
        resampling[[as.character(cnt)]] = tile
    }
    
    if (!return.factor) resampling = as.represampling(resampling)

    return(resampling)
}




# Function 'partition.kmeans'
# ---------------------------
# Uses k-means clustering to divide the samples
# into 'nfold' spatial clusters and perform
# spatial cross-validation.

#' Partition samples spatially using k-means clustering of the coordinates
#'
#' \code{partition.kmeans} divides the study area into irregularly shaped spatial partitions based on \emph{k}-means (\code{\link{kmeans}}) clustering of spatial coordinates.
#' @inheritParams partition.cv
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param nfold number of cross-validation folds, i.e. parameter \emph{k} in \emph{k}-means clustering
#' @param balancing.steps if \code{>1}, perform \code{nfold}-means clustering \code{balancing.steps} times, and pick the clustering that minimizes the Gini index of the sample size distribution among the partitions. The idea is that 'degenerate' partitions will be avoided, but this also has the side effect of reducing variation among partitioning repetitions. More meaningful constraints (e.g., minimum number of positive and negative samples within each partition should be added in the future.
#' @param order.clusters if \code{TRUE}, clusters are ordered by increasing x coordinate of center point
#' @param ... additional arguments to \code{\link{kmeans}} 
#' @return A \code{\link{represampling}} object, see also \code{\link{partition.cv}} for details.
#' @note Default parameter settings may change in future releases.
#' @references Brenning, A., S. Long & P. Fieguth. Forthcoming. Detecting rock glacier flow structures using Gabor filters and IKONOS imagery. Submitted to Remote Sensing of Environment.
#'
#' Russ, G. & A. Brenning. 2010a. Data mining in precision agriculture: Management of spatial information. In 13th International Conference on Information Processing and Management of Uncertainty, IPMU 2010; Dortmund; 28 June - 2 July 2010.  Lecture Notes in Computer Science, 6178 LNAI: 350-359.
#' @seealso \code{\link{sperrorest}}, \code{\link{partition.cv}}, \code{\link{partition.disc}}, \code{\link{partition.tiles}}, \code{\link{kmeans}}
#' @examples
#' data(ecuador)
#' resamp = partition.kmeans(ecuador, nfold = 5, repetition = 1:2)
#' plot(resamp, ecuador)
#' @export
partition.kmeans = function(data, coords = c("x","y"), nfold = 10,
    repetition = 1, seed1 = NULL, return.factor = FALSE, balancing.steps = 1,
    order.clusters = TRUE, ...)
{
    if (any(names(list(...)) == "kfold")) {
        warning("argument 'kfold' has been renamed to 'nfold' in 'partition.kmeans'")
        nfold = list(...)$kfold
    }
    balancing.steps = max(1, balancing.steps)

    resampling = list()
    for (cnt in repetition) {
        if (!is.null(seed1)) set.seed(seed1 + cnt)
        kms = list()
        for (i in 1:balancing.steps)
            kms[[i]] = kmeans(data[,coords], centers = nfold, ...)
        kmgini = function(x) {
            p = x$size / sum(x$size)
            return( 1 - sum(p^2) )
        }
        km = kms[[ which.max(sapply(kms, kmgini)) ]]
        # To do: add more meaningful selection criteria such as minimum number of positives and negatives in each partition ???
        if (order.clusters) {
            o = rank( km$center[,1], ties.method = "first" )
            km$cluster = o[km$cluster]
        }
        # The clusters are the partitions:
        tile = factor(km$cluster)
        
        if (!return.factor) tile = as.resampling(tile)
        resampling[[as.character(cnt)]] = tile
    }
    if (!return.factor) resampling = as.represampling(resampling)
    
    return(resampling)
}


#' Leave-one-disc-out cross-validation and leave-one-out cross-validation
#'
#' \code{partition.disc} partitions the sample into training and tests set by selecting circular 
#' test areas (possibly surrounded by an exclusion buffer) and using the remaining samples as training samples
#' (leave-one-disc-out cross-validation). \code{partition.loo} creates training and test sets for
#' leave-one-out cross-validation with (optional) buffer.
#' @name partition.disc
#' @inheritParams partition.cv
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param radius radius of test area discs; performs leave-one-out resampling if radius <0
#' @param buffer radius of additional 'neutral area' around test area discs that is excluded from training and test sets; defaults to 0, i.e. all samples are either in the test area or in the training area.
#' @param ndisc Number of discs to be randomly selected; each disc constitutes a separate test set. Defaults to \code{nrow(data)}, i.e. one disc around each sample.
#' @param return.train If \code{FALSE}, returns only test sample; if \code{TRUE}, also the training area.
#' @param prob optional argument to \code{\link{sample}}
#' @param replace optional argument to \code{\link{sample}}: sampling with or without replacement?
#' @param repetition see \code{partition.cv}; however, see Note below: \code{repetition} should normally be \code{=1} in this function.
#' @param ... arguments to be passed to \code{partition.disc}
#' @return A \code{\link{represampling}} object. Contains \code{length(repetition)} \code{resampling} objects. Each of these contains \code{ndisc} lists with indices of test and (if \code{return.train=TRUE}) training sets.
#' @note Test area discs are centered at (random) samples, not at general random locations. Test area discs may (and likely will) overlap independently of the value of \code{replace}. \code{replace} only controls the replacement of the center point of discs when drawing center points from the samples. \code{radius<0} does leave-one-out resampling with an optional buffer. \code{radius=0} is similar except that samples with identical coordinates would fall within the test area disc.
#' @references Brenning, A. 2005. Spatial prediction models for landslide hazards: review, comparison and evaluation. Natural Hazards and Earth System Sciences, 5(6): 853-862.
#' @seealso \code{\link{sperrorest}}, \code{\link{partition.cv}}, \code{\link{partition.kmeans}}
#' @examples
#' data(ecuador)
#' parti = partition.disc(ecuador, radius=200, buffer=200, ndisc=5, repetition=1:2)
#' plot(parti,ecuador)
#' summary(parti)
#' # leave-one-out with buffer:
#' parti.loo = partition.loo(ecuador, buffer=200)
#' summary(parti)
#' @export
partition.disc = function(data, coords = c("x","y"), radius, buffer = NULL, 
    ndisc = nrow(data), seed1 = NULL, return.train = TRUE,
    prob = NULL, replace = FALSE, repetition = 1)
{  
    posbuf = buffer
    if (is.null(buffer)) {
        pospuf = 0
    } else stopifnot(buffer >= 0)
    
    if (replace == FALSE & ndisc > nrow(data))
        stop("partition.disc: ndisc must be >nrow(data) if replace=FALSE")

    resample = list()

    # Loop for repetitions:
    for (cnt in repetition) {    
        if (!is.null(seed1)) set.seed(seed1 + cnt)
        if (ndisc == nrow(data)) {
            index = c(1:nrow(data))
        } else {
            index = sample.int(nrow(data), size=ndisc, replace = replace, prob = prob)
        }
        
        res = list()
        for (i in index) {
            if (!is.null(buffer) | radius >= 0) {
                di = sqrt( (data[,coords[1]] - data[i,coords[1]])^2 
                         + (data[,coords[2]] - data[i,coords[2]])^2 )
            }
            train.sel = numeric()
            if (radius >= 0) {
                # leave-disc-out with buffer:
                test.sel = which( di <= radius )
                if (return.train)
                    train.sel = which( di > (radius + posbuf) )
            } else {
                # leave-one-out with buffer:
                test.sel  = i
                if (return.train) {
                    if (is.null(buffer)) {
                        train.sel = c(1:nrow(data))[-i]
                    } else train.sel = which( di > posbuf )
                }
            }
            if (return.train & (length(train.sel) == 0)) {
                warning("empty training set in 'partition.disc': 'buffer' and/or 'radius' too large?")
            }
            res[[as.character(i)]] = list(train = train.sel, test = test.sel)
        }
        resample[[as.character(cnt)]] = res
    }    
    repres = as.represampling(resample)
    
    return(repres)
}


#' @rdname partition.disc
#' @name partition.loo
#' @export
partition.loo = function(data, ndisc = nrow(data), replace = FALSE, ...)
{
    partition.disc(data = data, radius = -1, ndisc = ndisc, replace = replace, ...)
}


#' Non-spatial bootstrap resampling
#'
#' \code{represampling.bootstrap} draws a bootstrap random sample (with replacement) from \code{data}.
#' @inheritParams partition.cv
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param nboot Size of bootstrap sample
#' @param oob logical (default \code{FALSE}): if \code{TRUE}, use the out-of-bag sample as 
#' the test sample; if \code{FALSE}, draw a second bootstrap sample of size \code{nboot} independently 
#' to obtain a test sample
#' @return A \code{\link{represampling}} object. This is a (named) list containing \code{length(repetition)}
#' \code{\link{resampling}} objects. Each of these contains only one list with indices of \code{train}ing and
#' \code{test} samples. Indices are row indices for \code{data}.
#' @examples
#' data(ecuador)
#' # only 10 bootstrap repetitions, normally use >=100:
#' parti = represampling.bootstrap(ecuador, repetition = 10)
#' plot(parti, ecuador) # careful: overplotting occurs 
#' # because some samples are included in both the training and
#' # the test sample (possibly even multiple times)
#' @export
represampling.bootstrap = function(data, coords = c("x","y"), nboot = nrow(data),
    repetition = 1, seed1 = NULL, oob = FALSE)
{
    resample = list()
    for (cnt in repetition) {
        if (!is.null(seed1)) set.seed(seed1 + cnt)
        # Bootstrap training sample, drawn with replacement:
        train = sample(nrow(data), nboot, replace = TRUE)
        if (oob) {
            # test set = out of bag sample:
            test = c(1:nrow(data))[ !(c(1:nrow(data)) %in% train) ]
        } else {
            # test set = independently drawn bootstrap sample
            test  = sample(nrow(data), nboot, replace = TRUE)
        }
        resample[[as.character(cnt)]] = list( "1" = list(train = train, test = test) )
    }
    return(as.represampling(resample))
}

#' Bootstrap at an aggregated level
#'
#' \code{represampling.factor.bootstrap} resamples partitions defined by a factor variable. This can be used for non-overlapping block bootstraps and similar.
#' @inheritParams represampling.bootstrap
#' @param fac defines a grouping or partitioning of the samples in \code{data}; three possible types: 
#' (1) the name of a variable in \code{data} (coerced to factor if not already a factor variable); 
#' (2) a factor variable (or a vector that can be coerced to factor); 
#' (4) a list of factor variables (or vectors that can be coerced to factor); this list must be of length \code{length(repetition)}, and if it is named, the names must be equal to \code{as.character(repetition)}; this list will typically be generated by a \code{partition.*} function with \code{return.factor=TRUE} (see Examples below)
#' @param nboot number of bootstrap replications used for generating the bootstrap training sample (\code{nboot[1]}) and the test sample (\code{nboot[2]}); \code{nboot[2]} is ignored (with a warning) if \code{oob=TRUE}. A value of \code{-1} will be substituted with the number of levels of the factor variable, corresponding to an \emph{n} out of \emph{n} bootstrap at the grouping level defined by \code{fac}.
#' @param oob if \code{TRUE}, the test sample will be the out-of-bag sample; if \code{FALSE} (default), the test sample is an independently drawn bootstrap sample of size \code{nboot[2]}
#' @details \code{nboot} refers to the number of groups (as defined by the factors) to be drawn with replacement from the set of groups. I.e., if \code{fac} is a factor variable, \code{nboot} would normally not be greater than \code{nlevels(fac)}, \code{nlevels(fac)} being the default as per \code{nboot=-1}.
#' @seealso \code{\link{represampling.disc.bootstrap}}, \code{\link{represampling.tile.bootstrap}}, note yet implemented: \code{partition.cv.factor}
#' @examples
#' data(ecuador)
#' # a dummy example for demonstration, performing bootstrap 
#' # at the level of an arbitrary factor variable:
#' parti = represampling.factor.bootstrap(ecuador, factor(floor(ecuador$dem/100)), oob=TRUE)
#' plot(parti,ecuador)
#' # using the factor bootstrap for a non-overlapping block bootstrap
#' # (see also represampling.tile.bootstrap):
#' fac = partition.tiles(ecuador, return.factor=TRUE, repetition=c(1:3), dsplit=500, min.n=200, rotation="random", offset="random")
#' parti = represampling.factor.bootstrap(ecuador, fac, oob=TRUE, repetition=c(1:3))
#' plot(parti,ecuador)
#' @export
represampling.factor.bootstrap = function(data, fac, repetition = 1, nboot = -1, seed1 = NULL, oob = FALSE)
{
    if (oob && length(nboot) > 1) warning("nboot[2] ignored because 'oob=TRUE'")
    if (is.list(fac)) {
        stopifnot(length(fac) == length(repetition))
        if (is.null(names(fac))) {
            names(fac) = as.character(repetition)
        } else stopifnot(all(as.character(repetition) %in% names(fac)))
    } else {
        if (length(fac) == 1 && is.character(fac)) {
            fac = data[,fac]
        } else stopifnot(length(fac) == nrow(data))
        fac = factor(fac)
    }
    if (length(nboot) == 1) nboot = rep(nboot,2)

    resample = list()

    for (cnt in repetition) {
        if (!is.null(seed1)) set.seed(seed1 + cnt)
        # what factor variable to resample?:
        if (is.list(fac)) {
            the.fac = factor(fac[[as.character(cnt)]])
        } else the.fac = fac
        # How many bootstrap samples (at the factor level)?:
        the.nboot = nboot
        if (the.nboot[1] < 0) the.nboot[1] = nlevels(the.fac)
        if (the.nboot[2] < 0) the.nboot[2] = nlevels(the.fac)
        # Factor levels to be used in training sample:
        train = sample(levels(the.fac), the.nboot[1], replace = TRUE)
        # Factor levels to be used for test sample:
        if (oob) { # out-of-bag, i.e. factors that are not used for the training sample:
            test = levels(the.fac)[ !(levels(the.fac) %in% train) ]
        } else {
            # second, independently drawn bootstrap sample at the factor level:
            test  = sample(levels(the.fac), the.nboot[2], replace = TRUE)
        }
        # Turn factor levels into indices:
        train = unlist( sapply(train, function(x) which(the.fac != x)), use.names = FALSE )
        test  = unlist( sapply(test,  function(x) which(the.fac == x)), use.names = FALSE )
        # Compile training and test indices into a resampling object:
        resample[[as.character(cnt)]] = as.resampling( list( "1" = list(train = train, test = test) ) )
    }
        
    resample = as.represampling(resample)
    return(resample)
}


#' Spatial block bootstrap using rectangular blocks
#'
#' \code{represampling.tile.bootstrap} performs a non-overlapping spatial block bootstrap by resampling at the level of rectangular partitions or 'tiles' generated by \code{\link{partition.tiles}}.
#' @inheritParams represampling.bootstrap
#' @param nboot see \code{\link{represampling.factor.bootstrap}}
#' @param oob see \code{\link{represampling.factor.bootstrap}}
#' @param ... additional arguments to be passed to \code{\link{partition.tiles}}
#' @export
represampling.tile.bootstrap = function(data, coords = c("x","y"), repetition = 1, nboot = -1, seed1 = NULL, oob = FALSE, ...)
{
    parti = partition.tiles(data = data, coords = coords, repetition = repetition, seed1 = seed1, return.factor = TRUE, ...)
    repres = represampling.factor.bootstrap(data = data, fac = parti, repetition = repetition, seed1 = seed1,
        nboot = nboot, oob = oob)
    return(repres)
}


#' Spatial block bootstrap at the level of spatial k-means clusters
#'
#' \code{represampling.kmeans.bootstrap} performs a non-overlapping spatial block bootstrap by resampling at the level of irregularly-shaped partitions generated by \code{\link{partition.kmeans}}.
#' @inheritParams represampling.bootstrap
#' @param nfold see \code{\link{partition.kmeans}}
#' @param nboot see \code{\link{represampling.factor.bootstrap}}
#' @param oob see \code{\link{represampling.factor.bootstrap}}
#' @param ... additional arguments to be passed to \code{\link{partition.kmeans}}
#' @export
represampling.kmeans.bootstrap = function(data, coords = c("x","y"), repetition = 1, nfold = 10, nboot = nfold, 
    seed1 = NULL, oob = FALSE, ...)
{
    parti = partition.tiles(data = data, coords = coords, repetition = repetition, seed1 = seed1, return.factor = TRUE, ...)
    repres = represampling.factor.bootstrap(data = data, fac = parti, repetition = repetition, seed1 = seed1,
        nboot = nboot, oob = oob)
    return(repres)
}


#' Overlapping spatial block bootstrap using circular blocks
#'
#' \code{represampling.disc.bootstrap} performs a spatial block bootstrap by resampling at the level of rectangular partitions or 'tiles' generated by \code{partition.tiles}.
#' @inheritParams represampling.bootstrap
#' @param oob logical (default \code{FALSE}): if \code{TRUE}, use the out-of-bag sample as
#' the test sample (the complement of the \code{nboot[1]} test set discs, minus the buffer area as specified in the \code{...} arguments to \code{\link{partition.disc}});
#' if \code{FALSE}, draw a second bootstrap sample of size \code{nboot} independently 
#' to obtain a test sample (sets of overlapping discs drawn with replacement)
#' @param nboot number of bootstrap samples; you may specify different values for the training sample (\code{nboot[1]}) and for the test sample (\code{nboot[2]})
#' @param ... additional arguments to be passed to \code{\link{partition.disc}}; note that a \code{buffer} argument has not effect if \code{oob=FALSE}; see example below
#' @note Performs \code{nboot} out of \code{nrow(data)} resampling of circular discs. This is an \emph{overlapping} spatial block bootstrap where the blocks are circular.
#' @examples
#' data(ecuador)
#' # Overlapping disc bootstrap:
#' parti = represampling.disc.bootstrap(ecuador, radius=200, nboot=20, oob=FALSE)
#' plot(parti,ecuador)
#' # Note that a 'buffer' argument would make no difference because boostrap sets of discs are
#' # drawn independently for the training and test sample.
#' #
#' # Overlapping disc bootstrap for training sample, out-of-bag sample as test sample:
#' parti = represampling.disc.bootstrap(ecuador, radius=200, buffer=200, nboot=10, oob=TRUE)
#' plot(parti,ecuador)
#' @export
represampling.disc.bootstrap = function(data, coords = c("x","y"), nboot,
    repetition = 1, seed1 = NULL, oob = FALSE, ...)
{
    if (oob && length(nboot) > 1) warning("nboot[2] ignored because oob=TRUE")
    if (length(nboot) == 1)  nboot = c(nboot, nboot)

    if (oob) {
        resample = list()
        for (cnt in repetition) {
            if (!is.null(seed1)) set.seed(seed1 + cnt)
            train = partition.disc(data = data, coords = coords, repetition = c(1:nboot[1]), 
                replace = TRUE, ndisc = 1, seed1 = NULL, return.train = TRUE, ...)
            test = c(1:nrow(data))
            for (i in 1:nboot[1]) {
                test = test[ test %in% train[[i]][[1]]$train ] # yes, $train!
                train[[i]][[1]]$train = NULL
            }
            if (length(test) == 0)
                warning("empty test set in 'partition.disc.bootstrap':\n'buffer' and/or 'radius' and/or 'nboot' too large?")
            train = unname(unlist(train))
            resample[[as.character(cnt)]] = list("1" = list(train = train, test = test))
        }
    } else {
        resample = list()
        for (cnt in repetition) {
            if (!is.null(seed1)) set.seed(seed1 + cnt)
            train = partition.disc(data = data, coords = coords, repetition = c(1:nboot[1]), seed1 = NULL, replace = TRUE, ndisc = 1, return.train = FALSE, ...)
            train = unname(unlist(train))
            test = partition.disc(data = data, coords = coords, repetition = c(1:nboot[2]), seed1 = NULL, replace = TRUE, ndisc = 1, return.train = FALSE, ...)
            test = unname(unlist(test))
            resample[[as.character(cnt)]] = list("1" = list(train = train, test = test))
        }
    }
        
    resample = as.represampling(resample)
    return(resample)
}



#' Plot spatial resampling objects
#'
#' \code{plot.represampling} displays the partitions or samples corresponding arising from the resampling of a data set
#' @method plot represampling
#' @name plot.represampling
#' @param x a \code{\link{represampling}} resp. \code{\link{resampling}} object
#' @param data a \code{data.frame} of samples containing at least the x and y coordinates of samples as specified by \code{coords}
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param pch point symbold (to be passed to \code{\link{points}})
#' @param wiggle.sd 'wiggle' the point locations in x and y direction to avoid overplotting of samples drawn multiple times by bootstrap methods;
#'  this is a standard deviation (in the units of the x/y coordinates) of a normal distribution and defaults to 0 (no wiggling)
#' @param ... additional arguments to \code{\link{plot}}
#' @note This function is not intended for samples obtained by resampling with replacement (e.g., bootstrap) 
#' because training and test points will be overplotted in that case. The size of the plotting region will also
#' limit the number of maps that can be displayed at once, i.e., the number of rows (repetitions) and
#' fields (columns).
#' @examples
#' data(ecuador)
#' # non-spatial cross-validation:
#' resamp = partition.cv(ecuador, nfold = 5, repetition = 1:2)
#' plot(resamp, ecuador)
#' # spatial cross-validation using k-means clustering:
#' resamp = partition.kmeans(ecuador, nfold = 5, repetition = 1:2)
#' plot(resamp, ecuador)
#' @export
plot.represampling = function(x, data, coords = c("x","y"), pch = "+", wiggle.sd = 0, ...)
{
    if (missing(data)) stop("'data' argument missing")
    stopifnot(wiggle.sd >= 0)
    resample = x
    nr = length(resample)
    nc = max( sapply(resample,length) )
    if (nr > 5) warning("Probably too many repetitions length(x) to be able to\nplot represampling object x. Trying anyway...")
    if (nr > 7) warning("Probably too many folds max(sapply(x,length)) to\nbe able to plot represampling object x. Trying anyway...")
    op = par(no.readonly = TRUE)
    par(mfrow = c(nr,nc), mar = c(2,2,3,.5), mgp = c(2,0.7,0), tcl = -0.3, cex = 0.5)
    for (i in 1:length(resample)) {
        for (j in 1:length(resample[[i]])) {
            seltrain = resample[[i]][[j]]$train
            seltest  = resample[[i]][[j]]$test
            main = paste("Repetition ", names(resample)[i], ", Fold ", j)
            plot(data[,coords[1]], data[,coords[2]], pch=".", type = "n",
                main = main, xlab = "", ylab = "", ...) #xlab=coords[1], ylab=coords[2])
            wxtrain = rnorm(length(seltrain), sd=wiggle.sd)
            wytrain = rnorm(length(seltrain), sd=wiggle.sd)
            wxtest = rnorm(length(seltest), sd=wiggle.sd)
            wytest = rnorm(length(seltest), sd=wiggle.sd)
            points(data[seltrain,coords[1]]+wxtrain, data[seltrain,coords[2]]+wytrain, pch=pch)
            points(data[seltest,coords[1]]+wxtest, data[seltest,coords[2]]+wxtest, pch=pch, col="red")
        }
    }
    par(op)
}

#' @rdname plot.represampling
#' @name plot.resampling
#' @method plot resampling
plot.resampling = function(x, ...)
{
    x = as.represampling(list( "1" = x ))
    plot.represampling(x)
}
