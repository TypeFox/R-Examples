#' @title Iterate over (N+1)D field and get all LC configurations
#' 
#' @description 
#' \code{data2LCs} gets all PLC or FLC configuration from a \eqn{(N+1)D} field 
#' given the LC template.  The shape and dimension of this LC template 
#' depends on coordinates passed on by \code{\link{setup_LC_geometry}}.
#' 
#' \subsection{User-defined LC template}{
#' 
#' Since \code{data2LCs} passes the \code{LC.coordinates} array to 
#' \code{\link{get_LC_config}} to iterate over the entire dataset, this 
#' functional programming approach allows user-defined light cone shapes 
#' (independent of the shapes implemented by \code{\link{setup_LC_geometry}}). 
#' 
#' Just replace the \code{$coordinates} from the \code{"LC"} class with a
#' user-specified LC template.
#' }
#' 
#' @param field spatio-temporal field; either a matrix or a 3-dimensional array 
#' with time \eqn{t} as the first dimension, and the spatial coordinates as 
#' subsequent dimensions. Make sure to check \code{\link{compute_LC_coordinates}} 
#' for correct formatting.
#' @param LC.coordinates coordinates for LC shape and dimension (usually the 
#' \code{$coordinates} value from the \code{"LC"} class; but also user-defined
#' coordinates are possible here).
#' @keywords manip
#' @export
#' @seealso \code{\link{compute_LC_coordinates}}, \code{\link{setup_LC_geometry}}
#' @examples
#' set.seed(1)
#' AA = matrix(rnorm(200), ncol = 10)
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 2, FLC = 0), shape ="cone")
#' bb = data2LCs(t(AA), LC.coordinates = LC_geom$coordinates)
#' image2(bb$PLC)
#' plot(density(bb$FLC))
#' 
#' # a time series example
#' data(nottem)
#' xx <- nottem
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 24, FLC = 3), space.dim = 0)
#' bb = data2LCs(xx, LC.coordinates = LC_geom$coordinates)
#' image2(bb$PLC)
#' plot(density(bb$FLC))
#' 
data2LCs <- function(field, 
                     LC.coordinates = list(PLC = NULL, FLC = NULL)){

  controls <- LC_coordinates2control_settings(LC.coordinates)
  n.p <- controls$n.p
  n.f <- controls$n.f
  space.dim <- controls$space.dim
  
  if (space.dim == 0) {
    if (length(dim(field)) > 0) {
      stop("The LC coordinates are for a time series; the field has at least one 
           space dimension. Please adjust.")
    }
  } else {
    if (space.dim + 1 != length(dim(field))){
      stop("The LC coordinates do not yield the same dimensionality for 
            the LC as the original field.")
    }
  }
  
  if (space.dim == 0) {
    #iter_spacetime <- get_spacetime_iterator(length(field), LC.coordinates)
    spacetime.grid <- get_spacetime_grid(length(field), LC.coordinates)
  } else {
    #iter_spacetime <- get_spacetime_iterator(dim(field), LC.coordinates)
    spacetime.grid <- get_spacetime_grid(dim(field), LC.coordinates)
  }
  LC_array <- matrix(NA, nrow = spacetime.grid$length.truncated,
                     ncol = 1 + 1 + space.dim + n.f + n.p )
  if (space.dim == 0) {
    colnames(LC_array) <- c("ID", "time", 
                           paste0("FLC", seq_len(n.f)), 
                           paste0("PLC", seq_len(n.p)))
  } else {
    colnames(LC_array) <- c("ID", "time", 
                            paste0("x", seq_len(space.dim)), 
                            paste0("FLC", seq_len(n.f)), 
                            paste0("PLC", seq_len(n.p)))
  }
  
  if (FALSE) {
    # works only with iterators; not available in R 3.0
    ii <- 0
    if (n.f == 1) {
      while (hasNext(iter_spacetime$truncated)){
        tmp.coords <- unlist(nextElem(iter_spacetime$truncated))
        ii <- ii + 1
        LC_array[ii, 1] <- ii
        LC_array[ii, 1 + seq_along(tmp.coords)] <- tmp.coords
        LC_array[ii, 1 + space.dim + 1 + seq_len(n.f)] <- c(field[tmp.coords[1], tmp.coords[2]])
        LC_array[ii, 1 + space.dim + 1 + n.f + seq_len(n.p)] <- get_LC_config(tmp.coords, field, LC.coordinates$PLC)
      }
    } else {
      while (hasNext(iter_spacetime$truncated)){
        tmp.coords <- unlist(nextElem(iter_spacetime$truncated))
        ii <- ii + 1
        LC_array[ii, 1] <- ii
        LC_array[ii, 1 + seq_along(tmp.coords)] <- tmp.coords
        LC_array[ii, 1 + space.dim + 1 + seq_len(n.f)] <- get_LC_config(tmp.coords, field, LC.coordinates$FLC)
        LC_array[ii, 1 + space.dim + 1 + n.f + seq_len(n.p)] <- get_LC_config(tmp.coords, field, LC.coordinates$PLC)
      }
    }
  }
  
  flc.cols <- 1 + 1 + space.dim + seq_len(n.f)
  plc.cols <- 1 + 1 + space.dim + n.f + seq_len(n.p)

  # work with tall matrix format grid values
  LC_array[, 1] <- seq_len(spacetime.grid$length.truncated)
  LC_array[, 1 + seq_len(space.dim + 1)] <- as.matrix(spacetime.grid$truncated)
  LC_array[, flc.cols] <-  
    t(get_LC_config(LC_array[, 1 + seq_len(space.dim + 1)], field, 
                    LC.coordinates$FLC))
  LC_array[, plc.cols] <- 
    t(get_LC_config(LC_array[, 1 + seq_len(space.dim + 1)], field, 
                    LC.coordinates$PLC))
  
  out <- list(FLC = cbind(LC_array[, flc.cols]),
              PLC = LC_array[, plc.cols])
 # out$dim <- list(original = iter_spacetime$dim.all, 
 #                  truncated = iter_spacetime$dim.truncated)
  out$dim <- list(original = spacetime.grid$dim.all, 
                  truncated = spacetime.grid$dim.truncated)
  if (length(dim(out$FLC)) == 0) {
    out$FLC <- matrix(out$FLC, ncol = 1)    
  }
  invisible(out)
}
