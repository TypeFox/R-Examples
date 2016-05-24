#' Creates the psychometric function
#' \code{create_psy_fun} creates the psychometric function
#' @keywords internal
#' @export
create_psy_fun <- function(psy_fun, guess, lapses) {
  psy_fun <- psy_fun
  function (x,p) {
    if (is.numeric(guess) && is.numeric(lapses)) {
      guess <- guess
      lapses <- lapses
      pshape <- p
    }
    if (is.logical(guess) && is.logical(lapses)){
      if (guess && lapses) {
        guess <- tail(p,2)[1]
        lapses <- tail(p,2)[2]
        pshape <- head(p,-2)
      }
      if (!guess && !lapses) {
        guess <- 0
        lapses <- 0
        pshape <- p
      }
    }
    if (is.numeric(guess) && is.logical(lapses)){
      if (lapses) {
        guess <- guess
        lapses <- tail(p,1)
        pshape <- head(p,-1)
      }
      if (!lapses) {
        guess <- guess
        lapses <- lapses
        pshape <- p
      }
    }
    if (is.logical(guess) && is.numeric(lapses)){
      if (guess) {
        guess <- tail(p,1)
        lapses <- lapses
        pshape <- head(p,-1)
      }
      if (!guess) {
        guess <- guess
        lapses <- lapses
        pshape <- p
      }
    }
    return(guess + (1 - guess- lapses) * psy_fun(x, pshape))
  }
}
