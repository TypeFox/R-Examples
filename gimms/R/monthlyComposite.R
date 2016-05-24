if ( !isGeneric("monthlyComposite") ) {
  setGeneric("monthlyComposite", function(x, ...)
    standardGeneric("monthlyComposite"))
}
#' Calculate monthly composite images
#'
#' @description
#' Based on a user-defined function, e.g. \code{max} for maximum value
#' composites (MVC), aggregate bi-monthly GIMMS datasets to monthly composites.
#'
#' @param x 'RasterStack' or 'character' vector of filenames. If the latter
#' applies and 'pos1', 'pos2' are not specified, the function will try to
#' retrieve monthly indices from \code{\link{monthlyIndices}}. Note that the
#' function does not work with binary data, but expects files that have
#' previously been created via \code{\link{rasterizeGimms}}.
#' @param indices 'numeric'. Indices identifying layers or files from identical
#' months.
#' @param fun Function. Used to calculate monthly composite layers, defaults to
#' \code{max}.
#' @param cores Integer. Number of cores for parallel computing.
#' @param filename Character. Optional output filename(s); see
#' \code{\link{writeRaster}}. If \code{cores > 1}, the number of supplied
#' filenames must match up with the number of unique monthly indices.
#' @param pos1,pos2 Numeric. If 'x' is a vector of filenames, the first and last
#' element of the date string to build monthly indices from. Defaults to the
#' GIMMS naming convention, see \code{\link{monthlyIndices}}.
#' @param ... Further arguments passed on to \code{\link{writeRaster}}.
#'
#' @return
#' A 'RasterStack' object with monthly composite layers.
#'
#' @author
#' Florian Detsch
#'
#' @seealso
#' \code{\link{stackApply}}, \code{\link{monthlyIndices}},
#' \code{\link{writeRaster}}.
#'
#' @examples
#' \dontrun{
#' ## Download sample data
#' gimms_dir <- paste0(getwd(), "/data")
#'
#' gimms_files <- downloadGimms(x = as.Date("2000-01-01"),
#'                              y = as.Date("2000-12-31"), dsn = gimms_dir)
#'
#' ## Rasterize files
#' gimms_raster <- rasterizeGimms(x = gimms_files, remove_header = TRUE)
#'
#' ## Create monthly maximum value composites
#' indices <- monthlyIndices(gimms_files)
#' gimms_raster_mvc <- monthlyComposite(gimms_raster, indices = indices)
#'
#' ## Visualize data
#' library(sp)
#' names(gimms_raster_mvc) <- paste(month.abb, 2000)
#' spplot(gimms_raster_mvc)
#' }
#'
#' @export monthlyComposite
#' @name monthlyComposite

################################################################################
### function using 'RasterStack' or 'RasterBrick' ##############################
#' @aliases monthlyComposite,RasterStackBrick-method
#' @rdname monthlyComposite
setMethod("monthlyComposite",
          signature(x = "RasterStackBrick"),
          function(x, indices, fun = max, cores = 1L, filename = "", ...) {

            ## stop if 'indices' is missing
            if (missing(indices))
              stop("Please supply a valid set of indices, e.g. returned by monthlyIndices().")

            ## check 'cores'
            cores <- checkCores(cores)

            ## immediately run 'stackApply'
            if (cores == 1L) {

              raster::stackApply(x, indices = indices, fun = fun,
                                 filename = filename, ...)

            ## or run it in parallel
            } else {

              # initialize cluster
              cl <- parallel::makePSOCKcluster(cores)
              doParallel::registerDoParallel(cl)

              # loop over unique layer indices and apply 'fun'
              i <- 1; j <- 1
              lst_out <- foreach::foreach(i = unique(indices)) %dopar% {
                x_sub <- raster::subset(x, which(indices == i))
                raster::stackApply(raster::subset(x, which(indices == i)),
                                   fun = fun,
                                   indices = rep(1, length(which(indices == i))),
                                   filename = ifelse(length(filename) == length(unique(indices)),
                                                     filename[which(unique(indices) == i)], ""),
                                   ...)
              }

              # deregister parallel backend
              parallel::stopCluster(cl)

              # stack layers
              raster::stack(lst_out)
            }
          })


################################################################################
### function using 'character' #################################################
#' @aliases monthlyComposite,character-method
#' @rdname monthlyComposite
setMethod("monthlyComposite",
          signature(x = "character"),
          function(x, pos1 = 4L, pos2 = 8L, fun = max, cores = 1L,
                   filename = "", ...) {

            ## check 'cores'
            cores <- checkCores(cores)

            ## extract timestamp from 'x'
            indices <- monthlyIndices(x, pos1 = pos1, pos2 = pos2)

            ## stack files and run 'monthlyComposite,RasterStackBrick-method'
            rst <- raster::stack(x)
            monthlyComposite(rst, indices = indices, fun = fun, cores = cores,
                             filename = filename, ...)

          })

