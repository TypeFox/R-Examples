#
# Copyright (c) 2010, Stephen B. Weston
#
# This is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

irep <- function(iterable, times, length.out, each) {
  # Apply "each" first
  it <- if (!missing(each)) {
    irep.each(iter(iterable), each)
  } else {
    iter(iterable)
  }

  if (!missing(length.out)) {
    # Ignore "times" if "length.out" is specified
    ilimit(recycle(it), length.out)
  } else if (!missing(times)) {
    if (length(times) == 1) {
      # If "times" has a single value, recycle that many times
      recycle(it, times)
    } else {
      # If "times" has multiple values, it's kind of like "each"
      irep.times(it, times)
    }
  } else {
    # Neither "length.out" or "times" was specified
    it
  }
}

# Internal function used to handle the irep "each" argument
irep.each <- function(it, each) {
  each <- as.integer(each[1])

  if (is.na(each)) {
    each <- 1L
  } else if (each < 0) {
    stop("invalid 'each' argument")
  }

  n <- 0L
  value <- NULL

  nextEl <- if (each == 0) {
    function() stop('StopIteration', call.=FALSE)
  } else if (each == 1) {
    function() nextElem(it)
  } else {
    function() {
      if (n <= 0) {
        value <<- nextElem(it)
        n <<- each
      }
      n <<- n - 1L
      value
    }
  }

  object <- list(nextElem=nextEl)
  class(object) <- c('abstractiter', 'iter')
  object
}

# Internal function used to handle the irep "times" argument
irep.times <- function(it, times) {
  times <- as.integer(times)
  if (length(times) == 0 || any(is.na(times) | times < 0)) {
    stop("invalid 'times' argument")
  }

  i <- 0L
  n <- 0L
  value <- NULL

  nextEl <- function() {
    while (n <= 0 && i < length(times)) {
      i <<- i + 1L
      n <<- times[i]
      value <<- nextElem(it)
    }
    if (n <= 0) {
      stop('StopIteration', call.=FALSE)
    }
    n <<- n - 1L
    value
  }

  object <- list(nextElem=nextEl)
  class(object) <- c('abstractiter', 'iter')
  object
}
