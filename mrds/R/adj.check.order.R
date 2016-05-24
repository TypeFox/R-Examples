#' Check order of adjustment terms
#'
#' 'adj.check.order' checks that the Cosine, Hermite or simple polynomials are
#' of the correct order.
#'
#'
#' Only even functions are allowed as adjustment terms. Also Hermite polynomials must be of degree at least 4 and Cosine of order at least 3. Finally, also checks that order of the terms >1 for half-normal/hazard-rate, as per p.47 of Buckland et al (2001). If incorrect terms are supplied then an error is throw via \code{stop}.
#'
#' @param adj.series Adjustment series used
#'   ('\code{cos}','\code{herm}','\code{poly}')
#' @param adj.order Integer to check
#' @param key key function to be used with this adjustment series
#' @return Nothing! Just calls \code{stop} if something goes wrong.
#'
#' @author David Miller
#' @seealso \code{\link{adjfct.cos}}, \code{\link{adjfct.poly}},
#'   \code{\link{adjfct.herm}}, \code{\link{detfct}}, \code{\link{mcds}},
#'   \code{\link{cds}}
#' @references S.T.Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake. 1993.
#'   Robust Models. In: Distance Sampling, eds. S.T.Buckland, D.R.Anderson,
#'   K.P. Burnham, J.L. Laake. Chapman & Hall.
#' @keywords methods
adj.check.order <- function(adj.series,adj.order,key){

  if(adj.series == "poly"){
    # If polynomial, check even in a very crude way
    if(any(as.integer(adj.order/2) != (adj.order/2))){
      stop("Odd polynomial adjustment terms selected")
    }

  }else if(adj.series == "herm"){
    # If hermite, check even and greater than (or equal to) order 4
    if(any(adj.order < 4)){
      stop("Hermite polynomial adjustment terms of order < 4 selected")
    }

    if(any(as.integer(adj.order/2) != (adj.order/2))){
      stop("Odd Hermite polynomial adjustment terms selected")
    }

  }

  adj.name <- switch(adj.series,
                     cos = "Cosine",
                     herm = "Hermite",
                     poly = "Simple polynomial",
                     NULL)

  key.name <- switch(key,
                     hn = "half-normal",
                     hr = "hazard-rate",
                     NULL)

  if((key %in% c("hn","hr")) & any(adj.order<2)){
    stop(paste0(adj.name," adjustments must be of order >2 for ",
                key.name," key functions"))
  }

  invisible()
}
