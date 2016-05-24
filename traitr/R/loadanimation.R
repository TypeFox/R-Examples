##  Copyright (C) 2010 John Verzani
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/

##' @include dialog.R
roxygen()

##' A window to show a loading animation
##'
##' @param message A message to display along with graphic while loading. PANGO markup is okay.
##' @return An item group instance with a \code{close} method to call to dismiss window
##' @export
##' @rdname misc
##' @examples
##' ## we call, something happens, then we close
##' \dontrun{
##' w <- loadingAnimation()
##' ## .... something long, like dlg$make_gui() ...
##' w$close()
##' }
loadingAnimation <- function(message="<b>Loading...</b>") {
  ## image from http://www.ajaxload.info/
  w <- gwindow("Loading", visible=FALSE, width=200, height=100)
  li <- anItemGroup(items=list(
                      labelItem(message, attr=c(markup=TRUE)),
                      imageItem(value=system.file("images/loading.gif", package="traitr"))
                      ),
                    w = w,
                    close = function(.) dispose(.$w)
                    )
  li$make_gui(container=w)
  visible(w) <- TRUE
  li
}
