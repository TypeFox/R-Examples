
#' Calculate mean nearest-neighbour distance between point datasets
#'
#' \code{dataset.distance} calculates Euclidean nearest-neighbour distances between two point datasets 
#' and summarizes these distances using some function, by default the mean.
#'
#' @param d1 a \code{data.frame} with (at least) columns with names given by \code{x.name} and \code{y.name}; these contain the x and y coordinates, respectively
#' @param d2 see \code{d1}  - second set of points
#' @param x.name name of column in \code{d1} and \code{d2} containing the x coordinates of points
#' @param y.name same for y coordinates
#' @param fun function to be applied to the vector of nearest-neighbor distances of \code{d1} from \code{d2}
#' @param method type of distance metric to be used; only \code{'euclidean'} is currently supported
#' @param ... additional arguments to \code{fun}
#' @return depends on \code{fun}; typically (e.g., \code{mean}) a numeric vector of length 1
#' @details Nearest-neighbour distances are calculated for each point in \code{d1}, resulting in a vector of length \code{nrow(d1)}, and \code{fun} is applied to this vector.
#' @examples
#' d = data.frame(x = rnorm(100), y = rnorm(100))
#' dataset.distance(d, d) # == 0
#' @export
#' @name dataset.distance
#' @aliases dataset.distance
#' @seealso \code{\link{add.distance}}
dataset.distance = function(d1, d2, 
    x.name="x", y.name="y", fun = mean, method = "euclidean", ...)
{
    method = tolower(method)
    if (method != "euclidean") {
        if (method == "euclidian") {
            warning("correct spelling is 'Euclidean', not 'Euclidian'")
            method = "euclidean"
        } else {
            warning("only Euclidean distance is currently implemented\n")
        }
    }
    di = rep(NA, nrow(d1))
    for (i in 1:nrow(d1))
        di[i] = min( sqrt( (d2[,x.name] - d1[i,x.name])^2 + (d2[,y.name] - d1[i,y.name])^2 ) )
    if (!is.null(fun)) di = fun(di,...)
    return( di )
}


#' Add distance information to resampling objects
#'
#' @name add.distance
#' @inheritParams partition.cv
#' @param object \code{\link{resampling}} or \code{\link{represampling}} object
#' @param ... Additional arguments to \code{\link{dataset.distance}} and \code{\link{add.distance.resampling}}, respectively
#' @return A \code{\link{resampling}} or \code{\link{represampling}} object containing an additional 
#' \code{$distance} component in each \code{\link{resampling}} object. The \code{distance} component is a
#' single numeric value indicating, for each \code{train} / \code{test} pair,
#' the (by default, mean) nearest-neighbour distance between the two sets.
#' @details Nearest-neighbour distances are calculated for each sample in the test set. These
#' \code{nrow(???$test)} nearest-neighbour distances are then averaged. Aggregation methods other
#' than \code{mean} can be chosen using the \code{fun} argument, which will be passed on to
#' \code{\link{dataset.distance}}.
#' @seealso \code{\link{dataset.distance}} \code{\link{represampling}} \code{\link{resampling}}
#' @examples
#' data(ecuador) # Muenchow et al. (2012), see ?ecuador
#' nsp.parti = partition.cv(ecuador)
#' sp.parti = partition.kmeans(ecuador)
#' nsp.parti = add.distance(nsp.parti, ecuador)
#' sp.parti = add.distance(sp.parti, ecuador)
#' # non-spatial partioning: very small test-training distance:
#' nsp.parti[[1]][[1]]$distance
#' # spatial partitioning: more substantial distance, depending on number of folds etc.
#' sp.parti[[1]][[1]]$distance
#' @export
add.distance = function(object, ...) UseMethod("add.distance")


#' @rdname add.distance
#' @name add.distance.resampling
#' @method add.distance resampling
#' @export
add.distance.resampling = function(object, data, coords = c("x", "y"), ...)
{
    for (j in 1:length(object)) {
        test.dist = dataset.distance(
            data[ object[[j]]$test, coords ], 
            data[ object[[j]]$train, coords ],
            x.name = coords[1], y.name = coords[2], ...)
        object[[j]]$distance = test.dist
    }
    return(object)
}

#' @rdname add.distance
#' @name add.distance.represampling
#' @method add.distance represampling
#' @export
add.distance.represampling = function(object, ...)
{
    object = lapply(object, add.distance.resampling, ...)
    class(object) = "represampling"
    return(object)
}

#' Alphanumeric tile names
#'
#' Functions for generating and handling alphanumeric tile names of the form \code{"X2:Y7"} as used by 
#' \code{\link{partition.tiles}} and \code{\link{represampling.tile.bootstrap}}.
#'
#' @name as.tilename
#' @aliases tilename
#' @param x object of class \code{tilename}, \code{character}, or \code{numeric} (of length 2)
#' @param ... additional arguments (currently ignored)
#' @return object of class \code{tilename}, \code{character}, or numeric vector of length 2
#' @examples
#' tnm = as.tilename(c(2,3))
#' tnm # "X2:Y3"
#' as.numeric(tnm) # c(2,3)
#' @seealso \code{\link{partition.tiles}}, \code{\link{represampling}}, \code{\link{represampling.tile.bootstrap}}
#' @export
as.tilename = function(x, ...) UseMethod("as.tilename")

#' @rdname as.tilename
#' @name as.tilename.numeric
#' @method as.tilename numeric
#' @export
as.tilename.numeric = function(x, ...) {
    stopifnot(length(x) == 2)
    stopifnot(is.numeric(x))
    stopifnot(all(x >= 0))
    stopifnot(x == round(x))
    x = paste("X", x[1], ":Y", x[2], sep="")
    class(x) = "tilename"
    return(x)
}

#' @rdname as.tilename
#' @name as.character.tilename
#' @method as.character tilename
#' @export
as.character.tilename = function(x, ...) {
    y = as.numeric.tilename(x) # just testing
    class(x) = "character"
    return(x)
}

#' @rdname as.tilename
#' @name as.numeric.tilename
#' @method as.numeric tilename
#' @export
as.numeric.tilename = function(x, ...) {
    x = strsplit(x,":")[[1]]
    if (length(x) != 2) stop("tilename objects must have the form 'X3:Y5' etc.")
    x = c( as.numeric(substring(x[1],2)), as.numeric(substring(x[2],2)) )
    stopifnot(all(!is.na(x)))
    stopifnot(all(x >= 0))
    return(x)
}

#' @rdname as.tilename
#' @name as.tilename.character
#' @method as.tilename character
#' @export
as.tilename.character = function(x, ...) {
    stopifnot(length(x) == 1)
    y = as.numeric.tilename(x) # just testing
    class(x) = "tilename"
    return(x)
}

#' @rdname as.tilename
#' @name print.tilename
#' @method print tilename
#' @export
print.tilename = function(x, ...) {
    print(as.character(x))
}




#' Identify small partitions that need to be fixed.
#'
#' \code{get.small.tiles} identifies partitions (tiles) that are too small according to some defined criterion / criteria (minimum number of samples in tile and/or minimum fraction of entire dataset).
#' @param tile factor: tile/partition names for all samples; names must be coercible to class \code{\link{tilename}}, i.e. of the form \code{"X4:Y2"} etc.
#' @param min.n integer (optional): minimum number of samples per partition
#' @param min.frac numeric >0, <1: minimum relative size of partition as percentage of sample
#' @param ignore character vector: names of tiles to be ignored, i.e. to be retained even if the inclusion criteria are not met.
#' @return character vector: names of tiles that are considered 'small' according to these criteria
#' @seealso \code{\link{partition.tiles}}, \code{\link{tilename}}
#' @examples
#' data(ecuador) # Muenchow et al. (2012), see ?ecuador
#' # Rectangular partitioning without removal of small tiles:
#' parti = partition.tiles(ecuador, nsplit = c(10,10), reassign = FALSE)
#' summary(parti)
#' length(parti[[1]])
#' # Same in factor format for the application of get.small.tiles:
#' parti.fac = partition.tiles(ecuador, nsplit = c(10,10), reassign = FALSE, return.factor = TRUE)
#' get.small.tiles(parti.fac[[1]], min.n = 20) # tiles with less than 20 samples
#' parti2 = partition.tiles(ecuador, nsplit = c(10,10), reassign = TRUE, 
#'     min.n = 20, min.frac = 0)
#' length(parti2[[1]]) # < length(parti[[1]])
#' @export
get.small.tiles = function(tile, min.n = NULL, min.frac = 0, ignore=c())
{
    stopifnot(is.factor(tile))
    # Number of samples in each tile:
    n.tile = tapply(tile, tile, length)
    # Number of tiles:
    n.tiles = length(n.tile)
    # Find the small ones:
    small.tile = rep(FALSE, n.tiles)
    if (is.null(min.n) & is.null(min.frac))
        stop("either 'min.n' or 'min.frac' must be specified in 'get.small.tiles'")
    if (!is.null(min.n))
        small.tile = small.tile | (n.tile < min.n)
    if (!is.null(min.frac))
        small.tile = small.tile | (n.tile*n.tiles/length(tile) < min.frac)
    if (any(small.tile)) {
        small.tile = levels(tile)[ small.tile ]
    } else small.tile = character()
    if ((length(small.tile) > 0) & (length(ignore) > 0))
        small.tile = small.tile[ !(small.tile %in% as.character(ignore)) ]
    # Order 'small' tiles, smallest one first:
    if (length(small.tile) > 0)
        small.tile = small.tile[ order(n.tile[small.tile], decreasing = FALSE) ]
    small.tile = factor(small.tile, levels = levels(tile))
    return(small.tile)
}

#' Determine the names of neighbouring tiles in a rectangular pattern
#'
#' This based on "counting" up and down based on the tile name.
#' @param nm Character string or factor: name of a tile, e.g., \code{"X4:Y6"}
#' @param tileset Admissible tile names; if missing and \code{nm} is a factor variable, then \code{levels(nm)} is used as a default for \code{tileset}
#' @param iterate internal - do not change default: to control behaviour in an interative call to this function
#' @param diagonal if \code{TRUE}, diagonal neighbours are also considered neighbours
#' @return Character string.
#' @name tile.neighbors
#' @export
tile.neighbors = function(nm, tileset, iterate = 0, diagonal = FALSE)
{
    if (missing(tileset)) {
        if (is.factor(nm)) {
            tileset = levels(nm)
        } else tileset = NULL
    }
    nm = as.character(nm)
    wh = as.numeric.tilename(as.tilename(nm))

    # Initial neighbors list:
    nbr = c()
    for (i in c(-1,0,1)) {
        if (wh[1] + i >= 0) {
            for (j in c(-1,0,1)) {
                if (!diagonal & i*j != 0) next
                if (i==j) next
                if (wh[2] + j >= 0) {
                    nbr = c( nbr, as.character( as.tilename( wh + c(i,j) ) ) )
                }
            }
        }
    }

    if (!is.null(tileset)) {
        # If neighbors are not in 'tileset' (i.e. have been removed e.g. because they were too small)
        # then continue search up to 'iterate' times/steps:
        if ( !any(nbr %in% tileset) & (iterate > 0) ) {
            nbrs = c()
            for (a.nbr in nbr)
                nbrs = c(nbrs, tile.neighbors(nm = a.nbr, tileset = tileset, iterate = iterate - 1))
            nbr = nbrs
            nbr = nbr[ nbr != nm ]
        }
        # Now remove any neighbor tiles that are not contained in the 'tileset', i.e. have
        # been removed previously (presumably because they were too small, or out of range):
        nbr = nbr[ nbr %in% tileset ]
    }

    if (!is.null(tileset))
        nbr = factor(nbr, levels = tileset)

    return(nbr)
}


#' Resampling objects such as partitionings or bootstrap samples
#'
#' Create/coerce and print resampling objects, e.g., partitionings or boostrap samples derived from a data set.
#'
#' @param object depending on the function/method, a list or a vector of type factor defining a partitioning 
#' of the dataset
#' @param x object of class \code{resampling}
#' @param ... currently not used
#' @name as.resampling
#' @return \code{as.resampling} methods: An object of class \code{resampling}.
#' @details A \code{resampling} object is a list of lists defining a set of training and test samples.
#'
#' In the case of \code{k}-fold cross-validation partitioning, for example, the corresponding \code{resampling} object
#' would be of length \code{k}, i.e. contain \code{k} lists. Each of these \code{k} lists defines a training
#' set of size \code{n(k-1)/k} (where \code{n} is the overall sample size), and a test set of size \code{n/k}.
#' The \code{resampling} object does, however, not contain the data itself, but only indices between \code{1} and \code{n}
#' identifying the seleciton (see Examples).
#'
#' Another example is bootstrap resampling. \code{\link{represampling.bootstrap}} with argument \code{oob=TRUE}
#' generates [\code{rep}]\code{resampling} objects with indices of a bootstrap sample in the \code{train} component
#' and indices of the out-of-bag sample in the test component (see Examples below).
#'
#' \code{as.resampling.factor}: For each factor level of the input variable, \code{as.resampling.factor} determines the indices of samples in this level (= test samples) and outside this level (= training samples). Empty levels of \code{object} are dropped without warning.
#'
#' \code{as.resampling.list} checks if the list in \code{object} has a valid \code{resampling} object structure (with components \code{train} and \code{test} etc.) and assigns the class attribute \code{"resampling"} if successful.
#' @examples
#' data(ecuador) # Muenchow et al. (2012), see ?ecuador
#'
#' # Partitioning by elevation classes in 200 m steps:
#' parti = factor( as.character( floor( ecuador$dem / 200 ) ) )
#' smp = as.resampling(parti)
#' summary(smp)
#' # Compare:
#' summary(parti)
#' 
#' # k-fold (non-spatial) cross-validation partitioning:
#' parti = partition.cv(ecuador)
#' parti = parti[[1]] # the first (and only) resampling object in parti
#' # data corresponding to the test sample of the first fold:
#' str( ecuador[ parti[[1]]$test , ] )
#' # the corresponding training sample - larger:
#' str( ecuador[ parti[[1]]$train , ] )
#'
#' # Bootstrap training sets, out-of-bag test sets:
#' parti = represampling.bootstrap(ecuador, oob = TRUE)
#' parti = parti[[1]] # the first (and only) resampling object in parti
#' # out-of-bag test sample: approx. one-third of nrow(ecuador):
#' str( ecuador[ parti[[1]]$test , ] )
#' # bootstrap training sample: same size as nrow(ecuador):
#' str( ecuador[ parti[[1]]$train , ] )
#' @seealso \code{\link{represampling}}, \code{\link{partition.cv}}, \code{\link{partition.kmeans}}, \code{\link{represampling.bootstrap}}, etc.
#' @aliases as.resampling resampling
#' @export
as.resampling = function(object, ...) {
    if (inherits(object,"resampling")) {
        return(object) 
    } else {
        return(UseMethod("as.resampling"))
    }
}


#' @rdname as.resampling
#' @name as.resampling.default
#' @method as.resampling default
#' @export
as.resampling.default = function(object, ...) {
    as.resampling.factor(factor(object))
}

#' @rdname as.resampling
#' @name as.resampling.factor
#' @method as.resampling factor
#' @export
as.resampling.factor = function(object, ...) {
    object = factor(object) # drop empty levels
    # Turn factor levels into test sets, one after the other:
    resampling = lapply(levels(object),
        function(x,spl) list( 
            train = which(spl != x), 
            test  = which(spl == x) ),
        spl = object)
    # result is a list with nlevels(object) levels
    names(resampling) = levels(object)
    class(resampling) = "resampling"
    return(resampling)
}

#' @rdname as.resampling
#' @name as.resampling.list
#' @method as.resampling list
#' @export
as.resampling.list = function(object, ...) {
    stopifnot(validate.resampling(object))
    class(object) = "resampling"
    return(object)
}

#' @rdname as.resampling
#' @name validate.resampling
#' @export
validate.resampling = function(object) {
    if (!is.list(object)) return(FALSE)
    for (i in 1:length(object)) {
        if (!is.list(object[[i]])) return(FALSE)
        if (length(object[[i]]) < 2) return(FALSE)
        if (!all(c("train","test") %in% names(object[[i]]))) return(FALSE)
        if (!is.numeric(object[[i]]$train) || !is.numeric(object[[i]]$test)) return(FALSE)
    }
    return(TRUE)
}

#' @rdname as.resampling
#' @name is.resampling
#' @export
is.resampling = function(x, ...) inherits(x, "resampling")

#' @rdname as.resampling
#' @name print.resampling
#' @method print resampling
#' @export
print.resampling = function(x, ...)
{
    cat("\nSample sizes in resampling object with", length(x), "folds:\n")
    print( as.data.frame( t( sapply( x, function(y) sapply(y, length) ) ) ) )
    cat("\n")
}


#' Resampling objects with repetition, i.e. sets of partitionings or boostrap samples
#'
#' Functions for handling \code{represampling} objects, i.e. \code{list}s of \code{\link{resampling}} objects.
#' @param object object of class \code{represampling}, or a list to be coerced to this class
#' @param x object of class \code{represampling}
#' @param ... currently not used
#' @name as.represampling
#' @return \code{as.represampling} methods return an object of class \code{represampling} with the 
#' contents of \code{object}.
#' @details \code{represampling} objects are (names) lists of \code{\link{resampling}} objects. Such
#' objects are typically created by \code{\link{partition.cv}}, \code{\link{partition.kmeans}},
#' \code{\link{represampling.disc.bootstrap}} and related functions.
#'
#' In \code{r}-repeated \code{k}-fold cross-validation, for example, the corresponding \code{represampling}
#' object has length \code{r}, and each of its \code{r} \code{\link{resampling}} objects has
#' length \code{k}.
#'
#' \code{as.resampling.list} coerces \code{object} to class \code{represampling} while coercing its elements 
#' to \code{\link{resampling}} objects. Some validity checks are performed.
#' @seealso \code{\link{resampling}}, \code{\link{partition.cv}}, \code{\link{partition.kmeans}}, \code{\link{represampling.disc.bootstrap}}, etc.
#' @examples
#' data(ecuador) # Muenchow et al. (2012), see ?ecuador
#' # Partitioning by elevation classes in 200 m steps:
#' fac = factor( as.character( floor( ecuador$dem / 300 ) ) )
#' summary(fac)
#' parti = as.resampling(fac)
#' # a list of lists specifying sets of training and test sets,
#' # using each factor at a time as the test set:
#' str(parti)
#' summary(parti)
#' @aliases as.represampling represampling
#' @export
as.represampling = function(object, ...)
    if (inherits(object,"represampling")) object else UseMethod("as.represampling")


#' @rdname as.represampling
#' @name as.represampling.list
#' @method as.represampling list
#' @export
as.represampling.list = function(object, ...) {
    valid = sapply(object, validate.resampling)
    if (any(!valid)) {
        msg = paste("cannot coerce to 'represampling' object: invalid list elements number\n   ", 
                paste(which(!valid),collapse=" ") )
        stop(msg)
    }
    object = lapply(object, as.resampling)
    class(object) = "represampling"
    return(object)
}

#' @rdname as.represampling
#' @name print.represampling
#' @method print represampling
#' @export
print.represampling = function(x, ...)
{
    txt = paste("Replicated Selection Object (r=", length(x), ")", sep = "")
    cat("\n", txt, "\n", sep = "")
    cat(paste(rep("-",nchar(txt)), collapse=""), "\n\n")
    for (i in 1:length(x)) print(x[[i]])
}

#' @rdname as.represampling
#' @name is.represampling
#' @export
is.represampling = function(object) inherits(object,"represampling")


#' Summary statistics for a resampling objects
#'
#' Calculates sample sizes of training and test sets within repetitions and folds of a \code{resampling} or \code{represampling} object.
#' @name summary.represampling
#' @method summary represampling
#' @param object A \code{resampling} or \code{represampling} object.
#' @param ... currently ignored
#' @return A list of \code{data.frame}s summarizing the sample sizes of training and test sets in each fold of each repetition
#' @export
summary.represampling = function(object, ...)
{
    lapply(object, function(x) as.data.frame( t(sapply(x, function(y) data.frame(n.train=length(y$train), n.test=length(y$test)))) ) )
}

#' @rdname summary.represampling
#' @name summary.resampling
#' @method summary resampling
#' @export
summary.resampling = function(object, ...)
{
    as.data.frame( t(sapply(object, function(y) data.frame(n.train=length(y$train), n.test=length(y$test)))) )
}
