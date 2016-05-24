#' Noise filtering through principal components
#' 
#' Filter noise from a RasterStack by decomposing into principal components 
#' and subsequent reconstruction using only a subset of components
#' 
#' @param x RasterStack to be filtered
#' @param k number of components to be kept for reconstruction 
#' (ignored if \code{expl.var} is supplied)
#' @param expl.var  minimum amount of variance to be kept after reconstruction
#' (should be set to NULL or omitted if \code{k} is supplied)
#' @param weighted logical. If \code{TRUE} the covariance matrix will be 
#' geographically weighted using the cosine of latitude during decomposition 
#' (only important for lat/lon data)
#' @param use.cpp logical. Determines whether to use \strong{Rcpp} 
#' functionality, defaults to \code{TRUE}.
#' @param verbose logical. If \code{TRUE} some details about the 
#' calculation process will be output to the console
#' @param ... additional arguments passed to \code{\link{princomp}}
#' 
#' @return a denoised RasterStack
#' 
#' @seealso
#' \code{\link{anomalize}}, \code{\link{deseason}}
#' 
#' @export denoise
#' 
#' @examples
#' data("vdendool")
#' vdd_dns <- denoise(vdendool, expl.var = 0.8)
#' 
#' opar <- par(mfrow = c(1,2))
#' plot(vdendool[[1]], main = "original")
#' plot(vdd_dns[[1]], main = "denoised")
#' par(opar)
denoise <- function(x,
                    k = NULL,
                    expl.var = NULL,
                    weighted = TRUE,
                    use.cpp = TRUE,
                    verbose = TRUE,
                    ...) {
  
  x.vals <- raster::getValues(x)
  #x.vals[is.na(x.vals)] <- 0
  
  # PCA
  if (weighted) { 
    pca <- princomp(~ x.vals, covmat = covWeight(x.vals, 
                                                 remote::getWeights(x)), 
                    scores = TRUE, na.action = na.exclude, ...)
  } else {
    pca <- princomp(~ x.vals, scores = TRUE, na.action = na.exclude, ...)
  }
  
  # declare reconstruction characteristics according to supplied values
  if (!is.null(expl.var)) {
    k <- which(cumsum(pca$sdev^2 / sum(pca$sdev^2)) >= expl.var)[1]
  } else {
    expl.var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))[k]
  }
  
  if (verbose) {
    cat("\n",
        "Using the first ",
        k,
        " components (of ",
        raster::nlayers(x),
        ") to reconstruct series...\n",
        " these account for ",
        expl.var,
        " of variance in orig. series\n\n", 
        sep = "")
  }
  
  # Reconstruction
  recons <- lapply(seq(nlayers(x)), function(i) {
    rowSums(t(as.matrix(pca$loadings[, 1:k])[i, ] * 
                t(pca$scores[, 1:k]))) + pca$center[i]
  })
  
  # Insert reconstructed values in original data set 
  if (use.cpp) { 
    mat <- raster::as.matrix(x)
    jnk <- insertReconsC(recons, mat)
    rst <- raster::setValues(x, jnk)
  } else {
    rst <- raster::brick(lapply(seq(recons), function(i) {
      tmp_rst <- x[[i]]
      tmp_rst[] <- recons[[i]]
      return(tmp_rst)
    }))
  }
  
  # Return denoised data set
  return(rst)
  
}
