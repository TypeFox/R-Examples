

#' Archetypes movies.
#' @param zs An \code{\link{archetypes}} object.
#' @param data The data matrix.
#' @param show Show archetypes or approximated data.
#' @param ssleep Seconds to sleep before start.
#' @param bsleep Seconds to sleep between each plot.
#' @param postfn Post plot function; is called in each
#'   iteration after the plot call.
#' @param rwdata.col1 If \code{show = 'rwdata'}: color of base data set.
#' @param rwdata.col2 If \code{show = 'rwdata'}: color of weighted data set.
#' @param ... Passed to underlying plot functions.
#' @return Undefined.
#' @aliases movieplot
#' @export
movieplot <- function(zs, data, show = c('atypes', 'adata', 'rwdata'),
                      ssleep = 0, bsleep = 0, postfn = function(iter){},
                      rwdata.col1 = gray(0.7), rwdata.col2 = 2, ...) {

  show <- match.arg(show)
  steps <- length(zs$history$states())

  if ( show == 'rwdata' )
    data <- zs$family$scalefn(t(data))

  # ... and play:
  Sys.sleep(ssleep)

  for ( i in seq_len(steps)-1 ) {
    a <- zs$history$get(i)[[1]]

    switch(show,
           atypes = {
             xyplot(a, data, ...)
           },
           adata = {
             plot(fitted(a), ...)
           },
           rwdata = {
             d <- zs$family$weightfn(data, a$reweights)

             plot(t(data), col = rwdata.col1, ...)
             points(t(d), col = rwdata.col2, ...)
           })

    postfn(i)

    Sys.sleep(bsleep)
  }
}


#' Archetypes plot movie 2.
#'
#' Shows the intermediate steps of the algorithm;
#'
#' @param zas.col Color of the intermediate archetypes.
#' @param zas.pch Type of the intermediate archetypes points.
#' @param old.col Color of the archetypes on step further.
#' @return Undefined.
#' @export
#' @rdname movieplot
movieplot2 <- function(zs, data, show='atypes',
                       ssleep=0, bsleep=0,
                       zas.col=2, zas.pch=13,
                       old.col=rgb(1,0.5,0.5), ...) {

  steps <- length(zs$history$states())

  Sys.sleep(ssleep)

  # Initial archetypes:
  a <- zs$history$get(0)[[1]]
  xyplot(a, data, ...)
  Sys.sleep(bsleep)

  # Alternating loop:
  for ( i in seq_len(steps-1) ) {
    a0 <- zs$history$get(i-1)[[1]]
    a <- zs$history$get(i)[[1]]

    xyplot(a0, data, atypes.col = old.col, ...)
    points(a$zas, col = zas.col, pch = zas.pch, ...)
    Sys.sleep(bsleep)

    xyplot(a0, data, atypes.col=old.col, ...)
    points(a$zas, col=zas.col, pch=zas.pch, ...)
    par(new=TRUE)
    xyplot(a, data, ...)
    Sys.sleep(bsleep)
  }

  xyplot(a, data, ...)
}


#' Archetypes parallel coordinates plot movie.
#' @return Undefined.
#' @export
#' @rdname movieplot
moviepcplot <- function(zs, data, show=c('atypes', 'adata'),
                        ssleep=0, bsleep=0, ...) {

  steps <- length(zs$history$states())
  atypesmovie <- ifelse(show[1] == 'atypes', TRUE, FALSE)
  rx <- apply(data, 2, range, na.rm=TRUE)

  Sys.sleep(ssleep)

  # ... and play:
  for ( i in seq_len(steps)-1 ) {
    a <- zs$history$get(i)

    if ( atypesmovie )
      pcplot(a, data, ...)
    else
      pcplot(fitted(a), rx=rx, ...)

    Sys.sleep(bsleep)
  }
}
