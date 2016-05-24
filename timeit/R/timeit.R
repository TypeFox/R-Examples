#' Profile a Function Call
#' 
#' This is a wrapper to \code{\link{Rprof}} that cleans up some of the
#' profile hand-holding and provides easier usage. This allows you to profile
#' either a single function call, or a whole block. Evaluation can be run
#' multiple times in order to assess variability in the timings of each
#' function call.
#' 
#' Function calls that get executed very quickly will be missed
#' by \code{Rprof}, unless you set \code{interval} very low. However, doing
#' this will probably break things (and isn't really important, since profiling
#' is there to help you catch the longest-running functions.) If you really
#' want to time quickly executed functions, you can manually set the \code{replications}
#' argument: we evaluate the \code{call} \code{replications} times, and infer
#' the (average) run-time of the function as \code{<time taken> / replications}.
#' 
#' @param call a call; this can be a single function call or a whole block.
#' @param replications integer; by default \code{NULL}, which indicates
#' we should 'guess' an appropriate number of replications. 
#' in order to more accurately profile
#' quickly-running functions, we run the call \code{replications} times,
#' and then infer the run-time as \code{<time>} / \code{replications}.
#' by default, the argument is \code{NULL} and we attempt to infer an
#' appropriate number of replications.
#' @param interval real. time interval between samples.
#' @param memory.profiling logical. output memory usage statistics?
#' @param times integer. how many times to call the function?
#' @param show.warnings boolean. output a warning if any iteration of the
#' run did not produce results?
#' @note If you set the \code{replications} argument high, you will likely
#' see some output from the \code{do_timeout} call that is unrelated to your
#' function call. This is due to all the wrapping of a function call to
#' be executed by \code{Rprof} introduces a minor overhead. For other
#' caveats, please see \code{\link{Rprof}}.
#' 
#' Currently, \code{timeit} does not support passing through of arguments,
#' so don't try to wrap \code{timeit} in a function call, whereby the 
#' call it attempts to evaluate is passed from a parent function. For example,
#' 
#' \code{f <- function(x) { timeit(x) }; f(rnorm(10))} 
#' 
#' won't work properly; a fix may come in the future.
#' 
#' @export
#' @return An object of S3 classes \code{timeit} and \code{data.frame}.
#' @seealso \code{\link{mean.timeit}} for mean running times over all
#' iterations processed, \code{\link{summary.timeit}} for summary
#' statistics,
#' \code{\link{plot.timeit}} for generating a boxplot of the returned
#' times, \code{\link{do_timeit}} for the workhorse function, and 
#' \code{\link{Rprof}} for information on how \R profiles the
#' execution of expressions.
#' @examples \dontrun{
#' tmp <- timeit({
#'   x <- 1:1E4; y <- x + runif(1E4)
#'   lm( y ~ x )
#'   }, times=5)
#' summary(tmp)
#' y <- 1E4
#' f <- function(x) { summary( sort( rnorm(x) ) ) }
#' tmp <- timeit( f(y), times=5 )
#' if( !is.null(tmp) ) {
#'   summary(tmp)
#'   mean(tmp)
#'   if( require(ggplot2) ) { plot(tmp) }
#' }}
timeit <- function(call,
                   replications=NULL,
                   interval=0.01,
                   memory.profiling=TRUE,
                   times=10,
                   show.warnings=FALSE
                   ) {
  
  ## arg pre-processing
  stopifnot( times > 0 )
  if( interval < 1E-3 ) {
    warning( "interval is very small; see ?Rprof for caveats on interval usage")
  }
  
  call_me <- match.call()$call
  if( length( grep( "^ ?call ?$", as.character(call_me), perl=TRUE ) ) > 0 ) {
    stop("'call' cannot be used as a variable / function name within your code block")
  }
  
  if( is.null(replications) ) {
    cat("Determining an appropriate number of replications... ")
    replications <- determine_replications( call_me, interval )
    cat("Done!\n", 
        replications, 
        if(replications==1) "replication" else "replications",
        "will be used.\n\n")
  }
  
  out_list <- vector("list", times)
  for( i in 1:times ) {
    cat("Running iteration", i, "of", times, "...\n")
    call_me <- match.call()$call
    out_list[[i]] <- do_timeit( call_me, 
                                replications=replications,
                                interval=interval, 
                                memory.profiling=memory.profiling,
                                show.warnings=show.warnings, 
                                i=i )
    
    ## re-extend out_list in case it has shrunk
    out_list[times*2] <- NULL
  }
  
  if( length( out_list ) == 0 ) {
    warning("No events were recorded. Try ",
            "setting 'replications' higher in the 'timeit' call.")
    return( invisible(NULL) )
  }
  
  out_list <- out_list[ sapply( out_list, function(x) { !is.null(x) } ) ]
  
  out <- as.data.frame( stringsAsFactors=FALSE, optional=TRUE,
                        do.call( rbind, out_list )
                        )
  out$func <- factor( out$func,
                      levels=names( sort( tapply( out$self.time, out$func, median ) ) )
                      )
  
  class(out) <- c("timeit", "data.frame")
  
  return(out)
  
}

#' Replication Helper
#' 
#' This function is used to infer an appropriate number of replications
#' in \code{\link{timeit}}.
#' 
#' @param call the call.
#' @param interval the timing interval as passed from \code{\link{timeit}}.
#' @param base the base number of replications to use.
determine_replications <- function( call, interval, base=1E6 ) {
  time <- microbenchmark( invisible( eval( call ) ), times=2 )$time[2]
  return( max( 1, as.integer( base / time / interval ) ) )
}

#' Profile a Function Call
#' 
#' This is the workhorse function called by \code{\link{timeit}}, and is
#' primarily meant to be called through \code{\link{timeit}}. However,
#' if you desire a more direct wrapper to \code{Rprof} then this can
#' be useful.
#' @param call a call (typically passed down through \code{timeit}).
#' @param interval real. time interval between samples.
#' @param replications integer; by default \code{NULL}, which indicates
#' we should 'guess' an appropriate number of replications. 
#' in order to more accurately profile
#' quickly-running functions, we run the call \code{replications} times,
#' and then infer the run-time as \code{<time>/replications}.
#' by default, the argument is \code{NULL} and we attempt to infer an
#' appropriate number of replications.
#' @param memory.profiling logical. include memory use in output?
#' @param show.warnings boolean. output a warning if any iteration of the
#' run did not produce results?
#' @param i integer. the iteration number. primarily for use from \code{\link{timeit}}.
#' @param gcFirst boolean. run the garbage collector before any evaluation of the function call?
#' @param gcDuring boolean. run the garbage collector before each iteration, as produced
#' by \code{replications}? (very slow)
#' @export
#' @return A data.frame of the profiling times.
do_timeit <- function(call, 
                      replications=NULL,
                      interval=0.005, 
                      memory.profiling=FALSE,
                      show.warnings=FALSE,
                      i=1,
                      gcFirst=TRUE,
                      gcDuring=FALSE
                      ) {
  
  if( !is.call(call) ) {
    call_me <- match.call()$call
  } else {
    call_me <- call
  }
  
  if( isTRUE( memory.profiling ) ) {
    memory <- "both"
  } else {
    memory <- "none"
  }
  
  if( is.null(replications) ) {
    replications <- determine_replications( call_me, interval )
  }
  
  if( gcFirst ) gc(FALSE)
  
  ..tmp.. <- tempfile()
  on.exit( unlink(..tmp..) )
  
  timeit_invisible <- invisible
  timeit_eval <- eval
  timeit_gc <- gc
  
  timeit_replicate <- function( replications, .call ) {
    for( i in 1:replications ) {
      if( gcDuring ) timeit_invisible( timeit_gc(FALSE) )
      timeit_invisible( timeit_eval( .call ) )
    }
  }
    
  Rprof( ..tmp.., interval=interval, memory.profiling=memory.profiling )
  timeit_replicate( replications, call_me )
  Rprof(NULL)
  
  out <- tryCatch( summaryRprof(..tmp.., memory=memory),
                   error = function(e) {
                     if( show.warnings ) warning("no events recorded for iteration ", i )
                     return( invisible(NULL) )
                   } )
  
  if( !is.null(out) ) {
    out <- out$by.self
    out$self.time <- out$self.time / replications
    out$total.time <- out$total.time / replications
    if( "mem.total" %in% names(out) ) {
      out$mem.total <- out$mem.total / replications
    }
    
    out$replications <- replications
    out$iteration <- i
    out$func <- rownames(out)
  
    ## remove the f'ns we don't need to know about
    if( length( grep( "timeit", rownames(out) ) ) > 0 ) {
      out <- out[ !(rownames(out) %in% grep( "timeit_", rownames(out), value=TRUE )), , drop=FALSE ]
    }
    
    if( nrow(out) == 0 ) {
      return( invisible(NULL) )
    } 
    
    ## re-calc the pct. times
    out$self.pct <- out$self.time / sum( out$self.time ) * 100
    out$total.pct <- out$total.time / sum( out$total.time ) * 100
    
  }
    
  return( out )
  
}

#' Summarize an 'timeit' Object
#' 
#' This function generates some summary statistics for a \code{\link{timeit}}
#' object.
#' 
#' @param object an object of class \code{timeit}.
#' @param ... unused additional arguments.
#' @export
#' @method summary timeit
#' @S3method summary timeit
summary.timeit <- function( object, ... ) {
  
  return(
    do.call( cbind, tapply( object$self.time, object$func, function(x) {
      tmp <- c( summary(x), length(x) )
      names(tmp)[7] <- "n"
      return(tmp)
    } ) )
  )
}

#' Plot a 'timeit' Object
#' 
#' This generates a boxplot of the timing output for a \code{timeit} object.
#' 
#' @param x the \code{timeit} object.
#' @param y unused.
#' @param min.pct number between 0 and 100. when set, we only plot
#' functions whose calling time makes up greater than \code{min.pct}
#' of the total calling time.
#' @param ... unused additional arguments.
#' @export
#' @method plot timeit
#' @S3method plot timeit
plot.timeit <- function( x, y=NULL, min.pct=5, ... ) {
  
  if( !require("ggplot2") ) {
    stop("plotting requires ggplot2")
  }
  
  if( min.pct < 0 || min.pct > 100 ) {
    stop( "min.pct must be between 0 and 100" )
  }
  
  num_iteration <- max( x$iteration )
  overall_times <- with( x, tapply( self.pct, func, function(x) {
    sum( x / num_iteration )
  } ) )
  keep <- names(overall_times)[ overall_times > min.pct ]
  x <- x[ x$func %in% keep, , drop=FALSE ]
  
  func <- as.symbol("func")
  self.time <- as.symbol("self.time")
  
  ## determine a time scale to use
  max.time <- max( x$self.time )
  if( max.time > 1 ) {
    use.time <- 1
    names(use.time) <- "seconds"
  } else if( max.time > 1E-3 ) {
    use.time <- 1E3
    names(use.time) <- "milliseconds"
  } else if( max.time > 1E-6 ) {
    use.time <- 1E6
    names(use.time) <- "microseconds"
  } else {
    use.time <- 1E9
    names(use.time) <- "nanoseconds"
  }
  
  x$self.time <- x$self.time * use.time
  
  print( ggplot( x, aes(x=func, y=self.time) ) +
           geom_boxplot(outlier.size=0) +
           geom_point( pch=21, fill="red", col="black", alpha=0.4 ) +
           xlab("") +
           ylab( paste( sep="", "Time (", names(use.time), ")" ) ) + 
           ggtitle( paste( sep="", "Time spent in each function call" ) ) +
           coord_flip()
  )
  
}

#' Calculate the mean of a 'timeit' Object
#' 
#' This function calculates the mean running time of each function call.
#' 
#' @param x the 'timeit' object.
#' @param ... additional arguments supplied to \code{\link{mean.default}}.
#' @export
#' @S3method mean timeit
#' @method mean timeit
mean.timeit <- function(x, ...) {
  names <- names(x)
  n <- max( x$iteration )
  fun <- function(x) {
    return( sum(x) / n )
  }
  out <- aggregate( x[ !(names(x) %in% "func") ], x["func"], FUN=fun )
  out$replications <- out$replications * n
  out$iterations <- tapply( x$iteration, x$func, length )
  out[ order(out$self.time, decreasing=TRUE), ]
}