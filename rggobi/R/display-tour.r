# Get tour projection
# Get the tour projection from a GGobi tour.
#
# This function retrieves the current projection matrix
# from a paused tour.  (The tour must be paused so that R
# can run commands).  
# 
# This can be used to record interesting projections of your
# data for later analysis.
#
# @arguments GGobiDisplay object running tour
# @keyword dynamic
#X g <- ggobi(mtcars)
#X d <- displays(g)[[1]]
#X \dontrun{
#X pmode(d) <- "2D Tour"
#X ggobi_display_get_tour_projection(d)
#X variables(d) <- list(X=names(mtcars))
#X ggobi_display_get_tour_projection(d)
#X MASS::eqscplot(as.matrix(mtcars) \%*\% ggobi_display_get_tour_projection(d))
#X }
ggobi_display_get_tour_projection <- function(gd) {
  mat <- .GGobiCall("getTourProjection", gd, pmode(gd))
  
  mat[,1:2] / mat[,3]
}

# Set tour projection
# Set the tour projection from a GGobi tour.
#
# If you know the projection you would like to see
# in the tour, you can use this function to set it.  The
# example illustrates setting the projection to show
# the first two principle components.
#
# @arguments GGobiDisplay object running tour
# @arguments tour projection
# @keyword dynamic
#X g <- ggobi(mtcars)
#X d <- displays(g)[[1]]
#X \dontrun{
#X pmode(d) <- "2D Tour"
#X variables(d) <- list(X=names(mtcars))
#X ggobi_display_get_tour_projection(d)
#X pc <- princomp(as.matrix(mtcars))$loadings[,1:2]
#X ggobi_display_set_tour_projection(d, pc)
#X pc <- princomp(as.matrix(mtcars), cor=T)$loadings
#X ggobi_display_set_tour_projection(d, pc)[,1:2]
#X }
ggobi_display_set_tour_projection <- function(gd, value) {
  normal <- all(abs(colSums(value^2) - 1) < 1e-3)
  orthog <- all(abs(crossprod(value, value) - diag(ncol(value))) < 1e-3)
  
  if (!normal) stop("Matrix is not normal (column lengths do not equal 1)")
  if (!orthog) stop("Matrix is not orthogonal")
  
  scale <- .GGobiCall("getTourProjection", gd, pmode(gd))[,3]
  value <- value * scale
  
  invisible(.GGobiCall("setTourProjection", gd, pmode(gd), value))
}