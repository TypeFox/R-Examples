#' Event Tables
#' 
#' Creates an event table, a custom \code{data.frame} used throughout the \code{linbin} package to store and manipulate linearly referenced data. Each row includes an event's endpoints \code{from} and \code{to} (which can be equal, to describe a point, or non-equal, to describe a line) and the values of any variables measured on that interval.
#' 
#' Event endpoints (and any additional arguments) are coerced to a data frame with \code{\link{data.frame}}, then coerced to an event table with \code{\link{as_events}}. A valid event table has two columns named "from" and "to" containing only finite numeric values (i.e., no \code{NA}, \code{NaN}, or \code{Inf}) and ordered such that \code{to} > or = \code{from}. \code{\link{is_events}} tests for these requirements. The other columns in the event table can be of any type supported by the \code{data.frame} class.
#' 
#' @param from,to Event endpoints, in any format coercible to single data frame columns. \code{from} and \code{to} are swapped as needed so that \code{to} > or = \code{from} on all rows. If \code{from} is the only non-empty argument, \code{\link{as_events}} is dispatched for object coercion.
#' @param ... Additional arguments, either of the form \code{value} or \code{tag = value}, to be passed directly to \code{\link{data.frame}} following \code{from} and \code{to}. Component names are created based on the tag (if present) or the deparsed argument itself.
#' @return An event table, the \code{data.frame} object used by \code{linbin} to describe interval data.
#' @seealso \code{\link{data.frame}}.
#' @seealso \code{\link{as_events}} and \code{\link{read_events}} for coercing objects and files to event tables, \code{\link{is_events}} to validate event tables.
#' @export
#' @examples
#' events(1, 5)
#' events(1:5)
#' events(c(0, 15, 25), c(10, 30, 35), x = 1, y = c('a', 'b', 'c'))
events <- function(from = numeric(), to = numeric(), ...) {
  # Attempt to coerce if only one non-empty argument
  if (!length(to) && !length(list(...)) && length(from)) {
    return(as_events(from))
  }
  # Otherwise, collect into data frame
  x <- data.frame(from, to, ...)
  return(as_events(x))
}

#' Coerce to an Event Table
#' 
#' Attempts to coerce an object to an event table.
#' 
#' @param x Object to be coerced to an event table.
#' @param from.col,to.col Names or indices of the columns in \code{x} containing the event endpoints. Values are swapped as needed to ensure that \code{to > or = from} on all rows.
#' @param ... Additional arguments passed to or used by methods.
#' @seealso \code{\link{events}} for creating event tables and \code{\link{read_events}} for reading files as event tables.
#' @export
#' @examples
#' as_events(1)
#' as_events(1:5)
#' as_events(cbind(1:5, 1:5), 1, 2)
#' as_events(data.frame(x = 1, start = 1:5, stop = 1:5), "start", "stop")
as_events <- function(x, ...) {
  if (is.null(x)) {
    return(data.frame(from = integer(), to = integer()))
  } else {
    UseMethod("as_events")
  }
}
#' @describeIn as_events Expands a numeric vector into two columns of event endpoints.
#' @export
as_events.numeric <- function(x, ...) {
  len <- length(x)
  if (len > 1) {
    # Interpret as sequential event endpoints
    e <- data.frame(from = x[-len], to = x[-1])
    need.flip <- e$from > e$to
    if (any(need.flip)) {
      e[need.flip, c("from", "to")] <- e[need.flip, c("to", "from")]
    }
    return(e)
  } else {
    # Interpret as single point event
    return(data.frame(from = x, to = x))
  }
}
#' @describeIn as_events Converts the matrix to a data frame, then calls the \code{data.frame} method.
#' @export
as_events.matrix <- function(x, from.col = 1, to.col = 2, ...) {
  return(as_events(as.data.frame(x), from.col = from.col, to.col = to.col, ...))
}
#' @describeIn as_events Renames \code{from.col} and \code{to.col} to "from" and "to" as needed. Since these column names must be unique, other columns cannot also be called "from" or "to".
#' @export
as_events.data.frame <- function(x, from.col = 1, to.col = 2, ...) {
  # Ensure endpoint columns exist and are unique
  if (is.character(from.col)) {
    from.col <- which(names(x) %in% from.col)
  }
  if (is.character(to.col)) {
    to.col <- which(names(x) %in% to.col)
  }
  names(x)[c(from.col[1], to.col[1])] <- c("from", "to")
  occurrence <- lapply(rgrep_exact(c("from", "to"), names(x)), length)
  if (any(occurrence < 1)) {
    stop("One or both of the required columns (from.col and to.col) are missing")
  }
  if (any(occurrence > 1)) {
    stop("One or both of the reserved column names ('from' and 'to') appear more than once")
  }
  # Coerce endpoint columns to numeric
  x[c("from", "to")] <- lapply(x[c("from", "to")], as.numeric)
  # Ensure endpoints are finite
  if (!all(is.finite(c(x$from, x$to)))) {
    stop("from.col and to.col cannot contain non-finite values (NA, NaN, and Inf)")
  }
  # Order endpoints (to > from)
  need.flip <- x$to < x$from
  if (any(need.flip)) {
    x[need.flip, c("from", "to")] <- x[need.flip, c("to", "from")]
  }
  return(x)
}

#' Read File as Event Table
#' 
#' Reads a file in table format and attempts to coerce it to an event table.
#' 
#' The file is read into R by calling \code{\link{read.table}}. Any of its arguments can be set by passing additional \code{tag = value} pairs. \code{from.col} and \code{to.col} are renamed to "from" and "to" as needed. Since these column names must be unique, other columns cannot also be called "from" or "to".
#' 
#' @param file Name, \code{\link{connection}}, or \code{\link{url}} of the file to be read as an event table.
#' @param from.col,to.col Names or indices of the columns containing event endpoints. Values are swapped as needed to ensure that \code{to > or = from} on all rows.
#' @param header Logical value indicating whether the file contains column names as its first line. If \code{FALSE}, columns will be named "V" followed by the column number, unless \code{col.names} (a vector of optional column names) is provided as an additional argument.
#' @param sep Character seperating values on each line of the file. If \code{sep = ""} (the default), the separator is 'white space' (that is, any combination of one or more spaces, tabs, newlines and carriage returns).
#' @param ... Additional arguments, of the form \code{tag = value}, to be passed directly to \code{\link{read.table}} to control how the file is read.
#' @seealso \code{\link{read.table}}.
#' @seealso \code{\link{events}} and \code{\link{as_events}} for creating event tables from existing objects.
#' @export
read_events <- function(file, from.col = 1, to.col = 2, sep = "", header = TRUE, ...) {
  x <- utils::read.table(file, sep = sep, header = header, ...)
  return(as_events(x, from.col, to.col))
}

#' Generate Sequential Events
#' 
#' Generates groups of regularly sequenced events fitted to the specified intervals. Intended for use as bins with \code{\link{sample_events}}.
#'
#' @param coverage An event table specifying the non-overlapping intervals to which the event sequences will be fitted. Gaps in coverage do not count towards event length. Points in the coverage are currently ignored.
#' @param length.out The number of events in each sequence. Event lengths are chosen such that they evenly divide the \code{coverage}.
#' @param by The length of the events in each sequence. Ignored if \code{length.out} is defined. When the length does not evenly divide the \code{coverage}, a shorter event is appended to the end of the sequence.
#' @param adaptive If \code{TRUE}, events are adjusted locally so that a whole number of events fit within each coverage interval, preserving breaks and gaps.
#' @return An endpoint-only event table with an additional group field if the length of \code{length.out} or \code{by} is \code{>} 1.
#' @seealso \code{\link{event_range}}, \code{\link{event_coverage}}, and \code{\link{fill_event_gaps}} for building a \code{coverage} from an existing event table.
#' @export
#' @examples
#' e <- events(c(0, 20, 40), c(10, 30, 45))
#' no.gaps <- event_range(e)
#' has.gaps <- event_coverage(e)
#' seq_events(no.gaps, by = 10)                           # unequal length (last is shorter)
#' seq_events(no.gaps, by = 10, adaptive = TRUE)          # equal length (11.25)
#' seq_events(no.gaps, length.out = 4)                    # equal length (11.25)
#' seq_events(has.gaps, length.out = 4, adaptive = FALSE) # equal coverage (11.25), straddling gaps
#' seq_events(has.gaps, length.out = 4, adaptive = TRUE)  # unequal coverage, fitted to gaps
#' seq_events(no.gaps, length.out = c(2, 4))              # "group" column added
seq_events <- function(coverage, length.out = NULL, by = NULL, adaptive = FALSE) {
  
  # Check inputs
  pts <- coverage$from == coverage$to
  if (any(pts)) {
    warning('ignoring points in coverage')
    coverage <- coverage[!pts, , drop = FALSE]
  }
  if (has_overlapping_events(coverage)) {
    stop('coverage cannot contain overlaps')
  }
  # Flatten event table, calculate total coverage
  total.length <- sum(coverage$to - coverage$from)
  if (!is.null(length.out))
    by <- total.length / round(as.numeric(length.out))
  
  # Generate bins for each bin length
  seq.bins <- lapply(by, function(by) {
    if (!adaptive) {
      # Build initial from and to values for the bins
      from <- min(coverage$from)
      to <- from + total.length
      binseq <- seq(from, to, by)
      if ((total.length / by) %% 1 != 0)
        # Add smaller bin to reach end of coverage
        binseq[length(binseq) + 1] <- to
      # Reinject gaps into bins
      gaps <- event_gaps(coverage)
      if (nrow(gaps)) {
        gap.length <- gaps$to - gaps$from
        # shift gaps by previous gaps' length
        gaps <- gaps - c(0, cumsum(gap.length)[-nrow(gaps)])
        # Locate start of gaps in bin sequence
        pos <- findInterval(gaps$from, binseq, rightmost.closed = TRUE)
        temp <- numeric(length(binseq))
        temp[unique(pos) + 1] <- stats::aggregate(gap.length, by = list(pos), sum)$x
        binseq <- binseq + cumsum(temp)
      }
      return(cbind(binseq[-length(binseq)], binseq[-1]))
    } else {
      # Fit bins to intervals of coverage
      seg.length <- coverage$to - coverage$from
      # Find evenly dividing length closest to nominal length
      r <- seg.length / by
      # Require at least one bin
      r[r < 1] <- 1
      l1 <- seg.length / ceiling(r)
      l2 <- seg.length / floor(r)
      d1 <- abs(l1 - by)
      d2 <- abs(l2 - by)
      smaller.d2 <- d2 < d1
      l1[smaller.d2] <- l2[smaller.d2]
      binmat <- do.call(rbind, lapply(seq_len(nrow(coverage)), function(i) {
        binseq <- seq(coverage$from[i], coverage$to[i], l1[i])
        return(cbind(binseq[-length(binseq)], binseq[-1]))
      }))
      return(binmat)
    }
  })
  # Format to event.table with group id
  if (length(seq.bins) > 1) {
    seq.bins <- Map(cbind, seq_along(seq.bins), seq.bins)
    seq.bins <- as_events(do.call(rbind, seq.bins), from.col = 2, to.col = 3)
    names(seq.bins)[1] = "group"
  } else {
    seq.bins <- as_events(do.call(rbind, seq.bins), 1, 2)
  }
  return(seq.bins)
}