#' @include archetypes-map.R
{}


#' Simplex visualization
#' 
#' The stochastic nature of the alpha coefficients implies that they 
#' exist on a standard (K-1)-simplex with the K archetypes Z as the
#' corners, and the coefficients as the coordinate with respect to these
#' corners. A standard simplex can be projected to two dimensions via 
#' a skew orthogonal projection, where all the vertices of the simplex
#' are shown on a circle connected by edges. The individual alpha 
#' coefficients can be then projected into this circle.
#' 
#' @param object An \code{\link{archetypes}} object
#' @param radius Radius of the projection
#' @param order Order of the archetypes
#' @param labels_cex Label expansion
#' @param labels Labels
#' @param show_labels Show labels
#' @param points_col Color of the points
#' @param points_pch Plot character of the points
#' @param points_cex Character expansion of the points
#' @param projection Projection function; see 
#'   \code{\link{archmap_projections}}
#' @param show_points Show the points
#' @param show_circle Show the circle
#' @param circle_col Color of the circle
#' @param show_edges Show the edges
#' @param edges_col Color of the edges
#' @param direction_length Expansion of the direction pointers
#' @param directions_col Color of the direction pointers
#' @param show_direction Show direction pointers
#' @param ... Additional arguments; currently ignored 
#'
#' @return
#'   Invisible list of all computed components needed for the simplex
#'   visualization.
#'   
#' @examples
#'   ### This example reproduces parts of the Figure 7 shown in 
#'   ### "Probabilistic Archetypal Analysis" by Seth and Eugster (2014)
#' 
#'   data("toy", package = "archetypes")
#'   
#'   set.seed(1234); a3 <- archetypes(toy, k = 3)
#'   set.seed(1237); a4 <- archetypes(toy, k = 4)
#'   set.seed(1238); a5 <- archetypes(toy, k = 5)
#'   
#'   simplexplot(a3)
#'   simplexplot(a3, show_direction = TRUE, show_points = FALSE)
#'   simplexplot(a4, projection = tspsimplex_projection)
#'   simplexplot(a5, show_direction = TRUE, show_points = FALSE, 
#'     direction_length = 2, directions_col = "black")
#'
#' @references
#'   See Section 6 in "Probabilistic Archetypal Analysis" by Seth and 
#'   Eugster (2014), http://arxiv.org/abs/1312.7604.   
#' 
#' @family simplexplot
#' 
#' @export
simplexplot <- function(object, radius = 10, order = NULL,
                        labels_cex = 1, labels = NULL, show_labels = TRUE,
                        points_col = "#00000044", points_pch = 19, points_cex = 1,
                        projection = simplex_projection, show_points = TRUE,
                        show_circle = TRUE, circle_col = "lightgray",
                        show_edges = TRUE, edges_col = "lightgray",
                        show_direction = FALSE,
                        direction_length = 1, directions_col = points_col, ...) {
  
  stopifnot("archetypes" %in% class(object))
  stopifnot(is.function(projection))
  
  k <- object$k
  
  if ( is.null(order) ) 
    order <- 1:k
  
  if ( is.null(labels) )
    labels <- sprintf("A%s", order)
  
  if ( length(points_col) == 1 )
    points_col <- rep(points_col, nrow(coef(object)))
  
  if ( length(points_cex) == 1 )
    points_cex <- rep(points_cex, nrow(coef(object)))
  
  if ( length(directions_col) == 1)
    directions_col <- rep(directions_col, nrow(coef(object)))
  
  
  params <- parameters(object)[order, ]
  coefs <- coef(object)[, order]
  
  
  proj_z <- projection(params, r = radius - 1)
  proj_h <- coefs %*% proj_z
  
  proj_labels <- proj_z
  t <- cbind(x = acos(proj_z[, "x"] / (radius-1)), y = asin(proj_z[, "y"] / (radius-1)))
  proj_labels <- cbind(x = radius * cos(t[, "x"]), y = radius * sin(t[, "y"]))
  
  proj_circle <- list(center = cbind(x = 0, y = 0), radius = radius - 1)
  proj_edges <- proj_z[as.integer(combn(1:k, 2)), ]
  
  proj_directions <- vector("list", length = nrow(object$alphas))

  for ( j in 1:nrow(object$alphas)) {
    s <- proj_h[j, , drop = FALSE]
    d <- matrix(NA_real_, ncol = 2, nrow = ncol(object$alphas))
    for ( i in 1:ncol(object$alphas) ) {  
      e <- proj_z[i, , drop = FALSE]
      
      v <- e - s
      m <- sqrt(sum(v^2))
      v <- v / m

      px <- s[1] + v[1] * direction_length * object$alphas[j, i]
      py <- s[2] + v[2] * direction_length * object$alphas[j, i]
      
      d[i, ] <- c(px, py)
    }
    proj_directions[[j]] <- list(s = s, e = d)
  }
  
  
  ### Plot:
  plot(proj_z, type = "n", asp = TRUE, 
       xlim = c(-radius, radius), ylim = c(-radius, radius),
       axes = FALSE, xlab = "", ylab = "")
  
  if ( show_circle ) {
    symbols(proj_circle$center, circles = radius - 1, 
            inches = FALSE, add = TRUE, asp = TRUE, fg = circle_col)
  }
  
  if ( show_edges ) {
    lines(proj_edges, col = edges_col)
  }
  
  if ( show_labels ) {
    text(proj_labels, labels = labels, cex = labels_cex)
  }
  
  if ( show_direction ) {
    for ( d in proj_directions ) {
      for ( i in 1:nrow(d$e) ) {
        lines(rbind(d$s, d$e[i, ,drop = FALSE]), col = directions_col) #[j])
      }
    }
  }
  
  if ( show_points ) {
    points(proj_h, col = points_col, pch = points_pch, cex = points_cex)
  }
  
  
  ret <- list(proj_z = proj_z, proj_h = proj_h, proj_labels = proj_labels,
              proj_directions = proj_directions, proj_circle = proj_circle,
              proj_edges = proj_edges)
  class(ret) <- "simplexplot"
 
  
  invisible(ret)
}



### Deviance: ########################################################

gaussian_deviance <- function(object, data) {
  y <- object$alphas %*% object$archetypes  
  sqrt(rowSums((y - data)^2))  
}


poission_deviance <- function(object, data) {
  t <- object$alphas %*% object$archetypes
  rowSums(2 * (data * log((data +.Machine$double.eps)/t) - data + t))
}


bernoulli_deviance <- function(object, data) {
  t <- object$alphas %*% object$archetypes
  t[t > 1] <- 1.0
  rowSums(2 * (data * log((data +.Machine$double.eps)/(t+.Machine$double.eps)) + 
                 (1-data) * log(((1-data) +.Machine$double.eps)/(1-t+.Machine$double.eps))))
}

