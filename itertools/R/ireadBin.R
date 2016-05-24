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

ireadBin <- function(con, what='raw', n=1L, size=NA_integer_,
                     signed=TRUE, endian=.Platform$endian, ipos=NULL) {
  # Sanity check "n"
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop('n must be a numeric value >= 1')
  }

  # Remember if we had to open this connection
  opened <- if (is.character(con)) {
    con <- file(con, open='rb')
    TRUE
  } else {
    if (!isOpen(con, 'r') || summary(con)$text != 'binary') {
      stop('con must be opened for reading in binary mode')
    }
    FALSE
  }

  if (!is.null(ipos)) {
    if (!isSeekable(con)) {
      stop('ipos cannot be specified unless con is seekable')
    }
    ipos <- iter(ipos)
  }

  nextEl <- function() {
    # Check if we've already stopped
    if (is.null(con)) {
      stop('StopIteration', call.=FALSE)
    }

    # "local" arguments to readBin that may be modified by "ipos"
    lwhat <- what
    ln <- n
    lsize <- size
    lsigned <- signed
    lendian <- endian

    # Seek on the connection if a position iterator has been specified
    if (!is.null(ipos)) {
      tryCatch({
        p <- nextElem(ipos)
      },
      error=function(e) {
        # Close the connection if necessary and propagate the exception
        if (opened) {
          close(con)
        }
        con <<- NULL
        stop(e)
      })

      # default value of "origin"
      origin <- 'start'

      if (is.list(p)) {
        # XXX should check for illegal element names in "p"
        # Don't do a "seek" unless a "where" value is specified
        if (!is.null(p$where)) {
          where <- p$where
          if (!is.null(p$origin)) origin <- p$origin
          seek(con, where=where, origin=origin, rw='read')
        }
        if (!is.null(p$what)) lwhat <- p$what
        if (!is.null(p$n)) ln <- p$n
        if (!is.null(p$size)) lsize <- p$size
        if (!is.null(p$signed)) lsigned <- p$signed
        if (!is.null(p$endian)) lendian <- p$endian
      } else {
        where <- p
        seek(con, where=where, origin=origin, rw='read')
      }
    }

    # Read the next "n" items
    d <- readBin(con, what=lwhat, n=ln, size=lsize, signed=lsigned, endian=lendian)

    # Check if we've hit EOF
    if (length(d) == 0) {
      # Close the connection if necessary
      if (opened) {
        close(con)
      }
      con <<- NULL
      stop('StopIteration', call.=FALSE)
    }

    d
  }

  it <- list(nextElem=nextEl)
  class(it) <- c('abstractiter', 'iter')
  it
}
