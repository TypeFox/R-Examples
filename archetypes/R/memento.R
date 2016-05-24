

#' Memento environment.
#'
#' Simple implementation of the 'Memento' design pattern.
#'
#' @param i The number of the state.
#' @param state The state to save.
#'
#' @return Memento environment.
#'
#' @examples
#'   \dontrun{
#'   m <- new.memento()
#'   m$save(i, state)
#'   m$states()
#'   m$get(i)
#'   }
#' @aliases memento
#' @rdname memento
#' @noRd
new.memento <- function() {

  memento <- new.env(parent = emptyenv())

  memento$save <- function(i, state) {
    assign(sprintf('s%s', i), state, envir = memento)
  }

  memento$get <- function(i) {
    if ( i < 0 )
      i <- length(memento$states()) + i - 1

    get(sprintf('s%s', i), envir = memento)
  }

  memento$states <- function() {
    ls(pattern = 's\\d+', envir = memento)
  }


  structure(memento, class = 'memento')
}
