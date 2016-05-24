#' Ripley's K-function
#' 
#' See \code{\link[spatstat]{Kest}} for more rigorous implementation in 2D.
#'
#' @param x Point pattern
#' @param r Vector of distances to estimate the function
#' @param ... ignored.
#' 
#' @details
#' Border correction is translation correction. 
#' Support in 3D therefore only for cuboidal windows.
#' 
#' @return
#' \code{\link{fv}}-object.
#' 
#' @useDynLib SGCS
#' @import spatstat
#' @export

Kfun <- function(x, r, ...) {
  ### prepare data
  x <- internalise_pp(x)
  r <- default_r(x, r)
  ### compute translation weights
  x$weights <- translation_weights(x)
  ### Compute:
  res <- .External("SGCS_Kfun_c",
                   x,
                   r,
                   PACKAGE="SGCS"
  )
  # scale:
  lambda <- x$n/x$area
  res <- res/lambda^2
  # make fv suitable
  C.final<-fv( data.frame(K=res, r=r, theo=rep(NA,length(r))),
               argu = "r",
               alim = range(r),
               ylab = substitute(K(r),NULL),
               desc = c("Ripley's K", "Poisson", "range"),
               valu = "K",
               fmla = ".~r",
               fname="K"
               )
  
  C.final  
}
