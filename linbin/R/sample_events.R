#' Find Intersecting Events
#' 
#' Returns a logical matrix indicating whether or not each pair of events intersect.
#' 
#' @param ex,ey Event tables.
#' @param equal.points If \code{TRUE}, equal-valued points are considered intersecting. This is always \code{TRUE} if \code{closed = TRUE}.
#' @param closed If \code{TRUE}, events are interpreted as closed intervals and events sharing only an endpoint are reported as intersecting.
#' @return A logical matrix with \code{ey} events as rows and \code{ex} events as columns.
#' @export
#' @examples
#' ex <- events(c(0, 5, 5, 10))
#' find_intersecting_events(ex, events(5), equal.points = FALSE) # equal points don't intersect
#' find_intersecting_events(ex, events(5), equal.points = TRUE)  # equal points do intersect
#' find_intersecting_events(ex, events(5), closed = TRUE)        # adjacent events intersect
#' find_intersecting_events(ex, ex)
find_intersecting_events <- function(ex, ey, equal.points = TRUE, closed = FALSE) {
  # FIXME: Slow, naive approach. Try speeding up by sorting first.
  inbins = apply(ex[c("from", "to")], 1, function(ex.each) {
    if (closed) {
      inbin <- ex.each[1] <= ey$to & ex.each[2] >= ey$from
    } else {
      inbin <- ex.each[1] < ey$to & ex.each[2] > ey$from
      if (equal.points) {
        inbin <- inbin | (ex.each[1] == ey$from & ex.each[2] == ey$to)
      }
    }
    return(inbin)
  })
  # Prevent returning a vector
  if (is.vector(inbins)) {
    inbins <- as.matrix(inbins)
  }
  return(inbins)
  # Converting from logical to a row column 
  # (for single matches): apply(inbins, 1, function(x) match(TRUE, x))
  # (for multiple matches): which(inbins, arr.ind = TRUE)
}

#' Sample Events
#' 
#' Computes event table variables over the specified sampling intervals, or "bins".
#' 
#' Events are cut at bin endpoints, and any \code{scaled.cols} columns are rescaled to the length of the resulting event segments. The event segments falling into each bin are passed to the sampling functions to compute the variables for each bin. Bins sample from events they overlap: line events with whom they share more than an endpoint, or point events with equal endpoints (if the bin itself is a point).
#' 
#' Sampling functions are specified in lists with the format \code{list(FUN, data.cols, by = group.cols, ...)}. The first element in the list is the function to use. It must compute a single value from one or more vectors of the same length. The following unnamed element is a vector specifying the event column names or indices to recursively pass as the first argument of the function. Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names. Additional unnamed elements are vectors specifying additional event columns to pass as the second, third, ... argument of the function. The first "by" element is a vector of event column names or indices used as grouping variables. Any additional named arguments are passed directly to the function. For example:
#' 
#' list(sum, 1:2, na.rm = TRUE) => sum(events[1], na.rm = TRUE), sum(events[2], na.rm = TRUE)
#' list(sum, 1, 3:4, 5) => sum(events[1], events[3], events[4], events[5]), ...
#' list(sum, c('x', 'y'), by = 3:4) => list(sum, 'x'), list(sum, 'y') grouped into all combinations of columns 3 and 4
#' 
#' Using the latter example above, column names are taken from the first argument (e.g. \code{x, y}), and all grouping variables are appended (e.g. \code{x.a, y.a, x.b, y.b}), where \code{a} and \code{b} are the levels of columns 3 and 4. \code{NA} is also treated as a factor level. Columns are added left to right in order of the sampling function arguments. Finally, names are made unique by appending sequence numbers to duplicates (using \code{\link{make.unique}}).
#' 
#' @param e An event table.
#' @param bins An event table specifying the intervals for sampling.
#' @param ... Lists specifying the sampling functions and parameters to be used (see the \code{Details}).
#' @param scaled.cols Names or indices of the event columns to be rescaled after cutting (see \code{\link{cut_events}}). Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names.
#' @param col.names Character vector of names for the columns output by the sampling functions. If \code{NULL}, the columns are named automatically (see the \code{Details}).
#' @param drop.empty If \code{TRUE}, bins not intersecting any events are dropped.
#' @return The \code{bins} event table with the columns output by the sampling functions appended.
#' @seealso \code{\link{seq_events}} to generate sequential bins.
#' @export
#' @examples
#' e <- events(from = c(0, 10, 15, 25), to = c(10, 20, 25, 40), length = c(10, 10, 10, 15), 
#'             x = c(1, 2, 1, 1), f = c('a', 'b', 'a', 'a'))
#' bins <- rbind(seq_events(event_coverage(e), 4), c(18, 18))
#' sample_events(e, bins, list(sum, 'length'))
#' sample_events(e, bins, list(sum, 'length'), scaled.cols = 'length')
#' sample_events(e, bins, list(sum, 'length', by = 'f'), scaled.cols = 'length')
#' sample_events(e, bins, list(weighted.mean, 'x', 'length'), scaled.cols = 'length')
#' sample_events(e, bins, list(paste0, 'f', collapse = "."))
sample_events <- function(e, bins, ..., scaled.cols = NULL, col.names = NULL, drop.empty = FALSE) {
  
  ### Initialize
  # Sampling functions
  functions <- sampling_functions(names(e), ...)
  # Calculate non-overlapping bin groups
  bin.groups <- group_nonoverlapping_events(bins)
  # Append original index to bins
  bins <- cbind(id = seq_len(nrow(bins)), bins)
  
  ### For each bin set
  L <- lapply(split(bins, bin.groups), function(bins) { 
    
    # Cut events at bin endpoints
    e.cut <- cut_events(e, bins, scaled.cols = scaled.cols)
    # Get index of segments intersecting with each bin
    inbins <- find_intersecting_events(bins, e.cut)
    # Grab first match for each event to convert to bin indices column
    # (each event in only one non-overlapping bin)
    # No match = NA
    bid <- apply(inbins, 1, function(x) match(TRUE, x))
    keep <- !is.na(bid)
    # Append assignments as last column, removing unassigned events
    e.assigned <- cbind(e.cut[keep, ], bin = bid[keep])
    # Apply the functions
    d <- do.call(cbind, lapply(functions, function(f) f(e.assigned)))
    # Reinsert empty bins as NA
    kept <- sort(unique(bid[keep]))
    if (drop.empty) {
      bins <- bins[kept, ]
    } else {
      nb <- nrow(bins)
      ind <- numeric(nb)
      ind[kept] <- seq_along(kept)
      ind[ind == 0] <- NA
      d <- d[ind, , drop = FALSE]
    }
    # Append bin info
    return(cbind(bins, d))
  })
  
  ### Prepare final output
  # Row bind list of outputs
  # FIXME: Can be slow for many/large output!
  # Convert to matrix (restoring factors) or use data.table (rbindlist)?
  # http://stackoverflow.com/questions/5980240/performance-of-rbind-data-frame
  D <- do.call(rbind, L)
  # Reorder result by original bin index and drop index
  # (has the effect of also making unique names)
  D <- D[order(D[[1]]), -1]
  # Apply user names
  if (!is.null(col.names)) {
    names(D)[(length(bins)):length(D)] <- col.names
  }
  return(D)
}

#' Build Function Call
#' 
#' Helper function for \code{\link{sampling_functions}}. Builds a function call in a enclosed environment with all fixed arguments and column indices so that the function can later be passed row subsets of an event table for sampling.
#' 
#' @param fun Function to use.
#' @param bin.col Column defining the groupping of bins.
#' @param data.cols Columns to each be passed as first argument to the function.
#' @param arg.cols Columns to be passed as the second, third, ... arguments of the function.
#' @param group.cols Columns to be used as factors.
#' @param arglist Additional arguments to pass to function \code{fun}.
#' @return A function whose parent environment encloses all of its fixed arguments.
#' @keywords internal
build_function_call <- function(fun, bin.col, data.cols, arg.cols = NULL, group.cols = NULL, arglist = NULL) {
  
  # Helpter function
  ident <- function(x) {
    y <- as.integer(as.factor(x))
    z <- gsub(" ", "0", format(y, scientific = FALSE))
    return(z)
  }
  
  ### Isolate function arguments from executing environment
  f <- new.env()
  f$fun <- fun
  f$bin.col <- bin.col
  f$data.cols <- data.cols
  f$arg.cols <- arg.cols
  f$group.cols <- group.cols
  f$arglist <- arglist
  
  ### Create function call where only input is the data
  f$call <- function(e) {
    
    # Build table of unique group combinations
    y <- e[c(bin.col, group.cols)]
    ynames <- names(y)
    if (ncol(y)) {
      grp <- rank(do.call(paste, c(lapply(rev(y), ident), list(sep = "."))), ties.method = "min")
    } else { 
      grp <- integer(NROW(e))
    }
    y <- y[match(sort(unique(grp)), grp, 0L), , drop = FALSE]
    
    # Pre-split additional argument columns
    var.args <- lapply(unname(e[arg.cols]), split, f = grp)
    
    # Apply function over each data column, then each subset of data and arguments
    e.data <- e[data.cols]
    z <- lapply(e.data, function(e.data.col) {
      args <- c(
        FUN = fun, 
        list((split(e.data.col, grp))),
        var.args,
        MoreArgs = list(arglist)
      )
      return(do.call(mapply, args))
    })
    
    # Assemble back into named dataframe
    len <- length(y)
    for (i in seq_along(z)) y[[len + i]] <- z[[i]]
    names(y) <- c(ynames, names(e.data))
    
    # Reshape long result to wide format
    # (if group cols present) 
    if (length(group.cols)) {
      # Merge group levels
      groups <- seq_along(group.cols) + 1
      y <- cbind(y, do.call(paste, c(y[groups], sep = ".")))
      # Reshape
      y <- stats::reshape(y[-groups], idvar = names(y[1]), timevar = names(y[length(y)]), direction = 'wide')
    }
    
    # Drop bin (id) column
    return(y[-1])
  }
  
  # Return function call
  return(f$call)
}

#' Build Sampling Functions
#' 
#' Helper function for \code{\link{sample_events}}. Parses function call parameters into self-enclosed function calls that can be passed row subsets of an event table for sampling.
#' 
#' NOTE: Assumes bin assignments will be appended to the end of the event table.
#' 
#' @param col.names Names of columns, for converting column name indices to numeric colum indices.
#' @param ... Lists of sampling function parameters (see \code{\link{sample_events}}).
#' @return A list of functions.
#' @keywords internal
sampling_functions <- function(col.names, ...) {
  
  ### Initialize inputs
  calls <- list(...)
  functions <- vector("list", length(calls))
  
  ### Parse each function call
  for (i in seq_along(calls)) {
    call <- calls[[i]]
    
    ## Pull out group.cols, if present
    group.ind <- match("by", names(call))
    if (!is.na(group.ind)) {
      group.cols <- rapply(call[group.ind], rgrep_exact, classes = c("character"), how = "replace", x = col.names)
      call[group.ind] <- NULL
    } else {
      group.cols <- NULL
    }
    args <- call[2:length(call)]
    ## Identify fixed arguments
    if (is.null(names(args))) 
      names(args) <- rep("", length(args))
    fixed <- !is.null(names(args)) & names(args) != ""
    
    ## Identify data and argument columns
    # Flatten to single level list (just in case)
    args[!fixed] <- lapply(args[!fixed], unlist)
    # Replace field names with corresponding numeric indices
    args[!fixed] <- lapply(rapply(args[!fixed], rgrep_exact, classes = "character", how = "replace", x = col.names), unlist)
    # Require integer numeric indices
    if (any(unlist(lapply(args[!fixed], is_not_integer))))
      stop('could not match names to column names')
    
    ## Build function call
    functions[[i]] <- build_function_call(
      fun = call[[1]], 
      bin.col = length(col.names) + 1,
      data.cols = args[!fixed][[1]],
      group.cols = unlist(group.cols),
      arg.cols = unlist(args[!fixed][-1]),
      arglist = args[fixed]
    )
  }
  
  ### Return list of all functions
  return(functions)
}