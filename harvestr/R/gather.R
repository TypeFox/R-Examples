{###############################################################################
# gather.R
# This file is part of the R package harvestr.
# 
# Copyright 2012 Andrew Redd
# Date: 6/2/2012
# 
# DESCRIPTION
# ===========
# Interface level functions that define the process flow.
#
# gather -> farm -> harvest
# gather -> plant -> harvest
# 
# LICENSE
# ========
# harvestr is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# dostats. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' Gather independent seeds.
#' 
#' Collect seeds representing independent random number streams.
#' These seeds can them be used in \code{\link{farm}} or \code{\link{plant}}.
#' 
#' 
#' @param x number of seeds, or an object with seeds to gather
#' @param seed a seed to use to set the seed, must be compatible with "L'Ecuyer-CMRG"
#' @param ... passed on
#' @param .starting if TRUE starting seeds with be gathered rather than ending seeds.
#' 
#' @seealso \link{RNG}
#' @family harvest
#' @importFrom parallel nextRNGStream
#' @include withseed.R
#' @export
gather <- 
function(x, seed=get.seed(), ..., .starting=F){
  oldseed <- get.seed()
  on.exit(replace.seed(oldseed))
  if(is.list(x)){
    seeds <- lapply(x, attr, ifelse(.starting,"starting.seed", "ending.seed"))
    if(any(sapply(seeds, is.null))) stop("Malformed list.")
    seeds
  } else if(is.numeric(x) && isTRUE(all.equal(x,ceiling(x)))){
    if(!is.null(seed)) {
        set.seed(seed, kind="L'Ecuyer-CMRG", normal.kind="Inversion")
    } else {
        RNGkind(kind="L'Ecuyer-CMRG", normal.kind="Inversion")
    }
    r <- get.seed()
    seeds <- vector('list', x)
    for(i in seq_len(x)) {
        r <-
        seeds[[i]] <-  structure(nextRNGStream(r), RNGlevel='stream')
    }
    structure(seeds, class=c('rng-seeds', 'list'))
  } else {
    stop("x must be either a list or integer")
  }
}

#' Create substreams of numbers based of a current stream.
#'
#' Seeds from \code{\link{gather}} can be used to generate another
#' set of independent streams.  These seeds can be given to 
#' \code{\link{graft}}
#' 
#' As a convenience \code{seed} can be an object that has a seed attached, 
#' ie. the result of any \code{harvestr} function.
#' 
#' @param seed a current random number stream compatible with 
#'        \code{\link{nextRNGSubStream}}
#' @param n number of new streams to create.
#'
#' @family harvest
#' @seealso \code{\link{nextRNGSubStream}}
#' @importFrom parallel nextRNGSubStream
#' @export
sprout <- function(seed, n) {
    es <- attr(seed, 'ending.seed')
    if(!is.null(es)) seed <- es
    rng.level <- attr(seed, 'RNGlevel')
    if(!is.null(rng.level) && rng.level == 'substream') {
        warning(paste('RNG seed provided si already a substream seed,'
                     ,'independence of streams not guaranteed.'))
    }
    
    seeds <- replicate(n, simplify=FALSE, {
        seed <<- structure(nextRNGSubStream(seed)
                          , RNGlevel = 'substream')
    }) 
    seeds
}

#' Call a function continuing the random number stream.
#' 
#' The \code{reap} function is the central function to \code{harvest}.
#' It takes an object, \code{x}, extracts the previous seed, ie. state of 
#' the random number generator, sets the seed, and continues any evaluation.
#' This creates a continuous random number stream, that is completly
#' reproducible.
#' 
#' The function calling works the same as the \link{apply} family of functions.
#' 
#' @param x         an object
#' @param fun       a function to call on object
#' @param ...       passed onto function
#' @param hash      hash of the list to retrieve the cache from.
#' @param cache     use cache, see Caching in \code{\link{harvestr}}
#' @param cache.dir directory for the cache.
#' @param time      should results be timed?
#' 
#' @seealso \code{\link{withseed}}, \code{\link{harvest}}, and \code{\link{with}}
#' @export
reap <-
function( x, fun, ...
        , hash      = digest(list(x, fun, ..., source="harvestr::reap"), "md5")
        , cache     = getOption('harvestr.use.cache', FALSE)
        , cache.dir = getOption("harvestr.cache.dir", "harvestr-cache")
        , time      = getOption('harvestr.time', FALSE)) {
  seed <- attr(x, "ending.seed")
  if(is.null(seed))
    stop("Could not find a seed value associated with x")
  if(cache){
    cache <- structure(cache, expr.md5 = hash)
  }
  f <- if(getOption('harvestr.use.try', TRUE)) {
    function(){try(fun(x,...), getOption('harvestr.try.silent', FALSE))}
  } else {
    function(){fun(x,...)}
  }
  withseed(seed, f, cache=cache, cache.dir=cache.dir, time=time)
}

#' Evaluate an expression for a set of seeds.
#' 
#' For each seed, set the seed, then evaluate the expression.
#' The \code{farm} function is used to generate data.
#' The Seeds for the state of the random numbrer generator is stored in the 
#' attribute \code{'ending.seed'}, and will be used by harvestr functions
#' for any other random number generation that is needed.
#' 
#' @param seeds a list of seeds can be obtained though \code{\link{gather}}
#' @param expr an expression to evalutate with the different seeds.
#' @param envir an environment within which to evaluate \code{expr}.
#' @param ... extra arguments
#' @param cache should cached results be used or generated?
#' @param time should results be timed?
#' @param .parallel should the computations be run in parallel?
#' 
#' @importFrom plyr llply ldply
#' @family harvest
#' @export
farm <-
function(seeds, expr, envir = parent.frame(), ...
        , cache     = getOption('harvestr.use.cache', FALSE)
        , time      = getOption('harvestr.time'     , FALSE)
        , .parallel = getOption('harvestr.parallel' , FALSE)){
  if(is.numeric(seeds) && length(seeds)==1)
    seeds <- gather(seeds)
  fun <- if(is.name(substitute(expr)) && is.function(expr)){
    stopifnot(is.null(formals(expr)))
    expr
  } else {
    substitute(expr)
  }
  results <- llply(.data=seeds, .fun=withseed, fun, envir=envir, ...
                  , cache=cache, time=time, .parallel=.parallel)
  if(time){
      attributes(results)$time <- total_time(results)
  }
  structure( results
           , 'function' = 'harvestr::farm'
           , class = c('harvestr-datasets', 'list')
           )
}


#' Harvest the results.
#' @param .list a list of \code{data.frame}s  See details.
#' @param fun a function to apply
#' @param ... passed to \code{fun}
#' @param time should results be timed?
#' @param .parallel should the computations be run in parallel?
#' 
#' @details
#' harvest is functionaly equivalant to llply, but takes on additional capability when used
#' with the other functions from this package.  When an object comes from \code{\link{withseed}}
#' the ending seed is extacted and used to continue evaluation.
#' 
#' @importFrom plyr mlply
#' @family harvest
#' @export
harvest <-
function(.list, fun, ...
        , time      = getOption('harvestr.time'     , FALSE)
        , .parallel = getOption('harvestr.parallel' , FALSE)){
  results <- llply(.list, reap, fun, ..., time =time,  .parallel=.parallel)
  if(getOption('harvestr.try.summary', TRUE)) try_summary(results)
  if(time){
      attributes(results)$time <- total_time(results)
  }
  structure( results
           , 'function' = 'harvestr::harvest'
           , class = c('harvestr-results', 'list')
           )
}

is_try_error <- function(x)inherits(x, 'try-error')
get_message <- function(e)attr(e, 'condition')$message
try_summary <- function(results){
    errors <- Filter(is_try_error, results)
    if(length(errors)){
        errors <- sapply(errors, get_message)
        cat(sprintf("%d of %d (%3.2g%%) calls produced errors"
                    , length(errors) , length(results)
                    , length(errors) / length(results) * 100), "\n\n")
        print(table(errors))
    }
}

#' Strip attributes from an object.
#' 
#' @param x, any object
#' @seealso \link{attributes}
#' @export
noattr <- noattributes <- function(x) {
  if(is.list(x)){
    x <- llply(x, noattributes)
  }
  attributes(x) <- NULL
  x
}

#' Assign elements of a list with seeds
#' @param .list a list to set seeds on
#' @param seeds to plant from \code{\link{gather}} or \code{\link{sprout}}
#'
#' @description
#' The function \code{plant} assigns 
#' each element in list set seed.
#' This will replace and ending seeds values already set for the 
#' objects in the list.
#' 
#' @family harvest
#' @export
plant <-
function(.list, seeds=gather(length(.list))) {
    if(inherits(.list, 'data.frame'))
        .list <- mlply(.list, data.frame)
    stopifnot(inherits(.list, 'list'))
    n <- length(.list)
    stopifnot(n == length(seeds))
    for(i in seq_len(n)){
        attr(.list[[i]],'ending.seed') <- seeds[[i]]
    }
    return(.list)
}

#' @rdname plant
#' @aliases graft
#' 
#' @description
#' The \code{graft} function replicates an object with independent substreams.
#' The result from \code{graft} should be used with \code{\link{harvest}}.
#' 
#' @param x an objects that already has seeds.
#' @param n number of seeds to create
#' @family harvest
#' @export
graft <-
function(x, n, seeds = sprout(x, n)) 
    plant(replicate(length(seeds), x, F), seeds)

    
#' Apply over rows of a data frame
#' 
#' @param df    a data frame of parameters
#' @param f     a function
#' @param ...   additional parameters
#' @param seeds seeds to use.
#' 
#' @return a list with f applied to each row of df.
#' 
#' @importFrom plyr splat
#' @export
plow  <-
function(df, f, ..., seeds=gather(nrow(df))){
    parameters <- plant(df, seeds=seeds)
    harvest(parameters, splat(f), ...)
}

#' Combine results into a data frame
#' 
#' @param l a list, from a harvestr function.
#' @param .check should checks be run on the object.
#' 
#' @seealso ldply
#' 
#' @export
bale <- 
function(l, .check=T){
    if(.check){
        stopifnot( is_homo(l)
                 , inherits(l, 'list')
                 , length(l)
                 )
        if(!inherits(l, 'harvestr::results'))
            warning('bale is intended to be used with harvestr results but got a ', class(l))
    }
    ldply(l, if(inherits(l[[1]], 'harvestr-results')) bale else I)
}




