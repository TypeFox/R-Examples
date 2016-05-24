#
# Copyright (c) 2012, Stephen B. Weston
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

irecord <- function(con, iterable) {
  if (is.character(con)) {
    con <- file(con, 'wb')
    on.exit(close(con))
  }
  it <- iter(iterable)
  tryCatch({
    repeat {
      serialize(nextElem(it), con)
    }
  },
  error=function(e) {
    if (! identical(conditionMessage(e), 'StopIteration'))
      stop(e)
  })
  invisible()
}

ireplay <- function(con) {
  # Remember if we had to open this connection
  opened <- if (is.character(con)) {
    con <- file(con, open='rb')
    TRUE
  } else {
    FALSE
  }

  nextEl <- function() {
    # Check if we've already stopped
    if (is.null(con)) {
      stop('StopIteration', call.=FALSE)
    }

    tryCatch({
      unserialize(con)
    },
    error=function(e) {
      if (opened) {
        close(con)
      }
      con <<- NULL
      stop('StopIteration', call.=FALSE)
    })
  }

  object <- list(nextElem=nextEl)
  class(object) <- c('abstractiter', 'iter')
  object
}
