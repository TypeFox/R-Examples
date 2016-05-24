#
# Copyright (c) 2009-2013, Stephen B. Weston
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

ichunk <- function(iterable, chunkSize, mode='list') {
  force(iterable)
  force(chunkSize)
  it <- iter(iterable)

  legal.modes <- c('list', 'logical', 'integer', 'numeric', 'double',
                   'complex', 'character', 'raw')
  if (! mode %in% legal.modes)
    stop(sprintf("cannot make a vector of mode '%s'", mode))

  nextEl.list <- function() {
    r <- vector('list', chunkSize)
    i <- 0L

    tryCatch({
      while (i < chunkSize) {
        r[i + 1L] <- list(nextElem(it))
        i <- i + 1L
      }
    },
    error=function(e) {
      if (!identical(conditionMessage(e), 'StopIteration') || i == 0L)
        stop(e)
      length(r) <<- i
    })

    r
  }

  nextEl.vector <- function() {
    r <- vector(mode, chunkSize)
    i <- 0L

    tryCatch({
      while (i < chunkSize) {
        r[i + 1L] <- nextElem(it)
        i <- i + 1L
      }
    },
    error=function(e) {
      if (!identical(conditionMessage(e), 'StopIteration') || i == 0L)
        stop(e)
      length(r) <<- i
    })

    r
  }

  object <- if (mode == 'list')
    list(nextElem=nextEl.list)
  else
    list(nextElem=nextEl.vector)
  class(object) <- c('abstractiter', 'iter')
  object
}
