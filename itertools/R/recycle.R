#
# Copyright (c) 2009-2010, Stephen B. Weston
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

recycle <- function(iterable, times=NA_integer_) {
  # Manually check for a missing argument since "inherits" issues
  # a cryptic error message in that case
  if (missing(iterable)) {
    stop('argument "iterable" is missing, with no default')
  }

  if (!is.numeric(times) || length(times) != 1 || (!is.na(times) && times < 0)) {
    stop('argument "times" must be a non-negative numeric value')
  }

  times <- as.integer(times)

  if (is.na(times) || times > 1) {
    if (! inherits(iterable, 'iter')) {
      buffer <- iterable
      buffer.iter <- iter(buffer)
    } else {
      iterable.iter <- iter(iterable)
      bsize <- 256  # allocated size of buffer
      bsize.max <- 2 ^ 31 - 1  # maximum allowable allocated size of buffer
      buffer <- vector('list', length=bsize)
      blen <- 0  # number of values currently in buffer
      buffer.iter <- NULL  # will become an iterator over buffer
    }
  } else if (times > 0) {
    iterable.iter <- iter(iterable)
  }

  # This is used until the underlying iterator runs out
  nextEl.buffering <- function() {
    tryCatch({
      # Check if buffer is full
      if (blen >= bsize) {
        # Don't attempt to create a list with more than 2^31-1 elements
        if (blen == bsize.max) {
          stop('underlying iterator has too many values to buffer')
        }
        # Double the size of buffer
        bsize <<- min(2 * bsize, bsize.max)
        length(buffer) <<- bsize
      }
      e <- nextElem(iterable.iter)
      blen <<- blen + 1
      buffer[blen] <<- list(e)
      e
    },
    error=function(e) {
      if (identical(conditionMessage(e), 'StopIteration')) {
        times <<- times - 1L  # will still be greater than zero
        length(buffer) <<- blen
        iterable <<- NULL
        iterable.iter <<- NULL
        buffer.iter <<- iter(buffer)
        nextEl.pointer <<- nextEl.cycling
        nextEl()
      } else {
        stop(e)
      }
    })
  }

  # This will be used once we've run through the underlying iterator
  nextEl.cycling <- function() {
    tryCatch({
      nextElem(buffer.iter)
    },
    error=function(e) {
      if (identical(conditionMessage(e), 'StopIteration')) {
        if (!is.na(times) && times <= 1) {
          times <<- 0L
          stop(e)
        }
        times <<- times - 1L
        buffer.iter <<- iter(buffer)
        # If this throws 'StopIteration', we're done
        nextElem(buffer.iter)
      } else {
        stop(e)
      }
    })
  }

  # This handles the case when "times" is one (pretty useless case)
  nextEl.one <- function() {
    nextElem(iterable.iter)
  }

  # This handles the case when "times" is zero
  nextEl.zero <- function() {
    stop('StopIteration', call.=FALSE)
  }

  # Set the initial value of nextEl.pointer
  if (is.na(times) || times > 1) {
    nextEl.pointer <- if (is.null(buffer.iter)) nextEl.buffering else nextEl.cycling
  } else if (times == 1) {
    nextEl.pointer <- nextEl.one
  } else {
    nextEl.pointer <- nextEl.zero
  }

  # This is the function that will be stored in the iterator object,
  # which will call either nextEl.buffering of nextEl.cycling, depending
  # on the value of nextEl.pointer variable
  nextEl <- function() {
    nextEl.pointer()
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}
