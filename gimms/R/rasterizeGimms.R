#' Rasterize GIMMS NDVI3g binary data
#'
#' @description
#' Import GIMMS NDVI3g binary data into R as 'Raster*' objects based on a
#' companion header file.
#'
#' @param x Character. Vector of local filepaths.
#' @param flag Logical. If \code{TRUE}, returns flag values instead of NDVI. See
#' 'References' for further reading. Note that this will be ignored if
#' \code{scaling = FALSE}.
#' @param header Character. Header file(s) corresponding to the binary data
#' in 'x'. If missing, the standard header file for GIMMS NDVI3g binary data
#' will be used.
#' @param water2na Logical. Determines whether or not to discard pixels with
#' 'mask-water' value (-10000; see 'References').
#' @param nodata2na Logical. Determines whether or not to discard pixels with
#' 'mask-nodata' value (-5000; see 'References').
#' @param scaling Logical. If \code{TRUE} (default), scaling is enabled which
#' allows to retrieve NDVI or flag values (depending on 'flag').
#' @param remove_header Logical. If \code{TRUE}, the header file specified in
#' 'header' or, if not specified, created internally via
#' \code{gimms:::createHeader} is removed after all operations have finished.
#' @param cores Integer. Number of cores for parallel computing. If 'filename'
#' is not specified, the parallel option is automatically disabled.
#' @param filename Character. Optional vector of output filenames with the same
#' length as 'x'; see \code{\link{writeRaster}}.
#' @param ... Further arguments passed on to \code{\link{writeRaster}} (except
#' for 'bylayer' and 'suffix').
#'
#' @return
#' If 'x' is a single filename, an object of class 'RasterLayer';
#' if 'x' is a vector of filenames, an object of class 'RasterStack'.
#'
#' @author
#' Florian Detsch
#'
#' @seealso
#' \code{gimms:::createHeader}, \code{\link{raster}}, \code{\link{writeRaster}}.
#'
#' @references
#' \url{https://nex.nasa.gov/nex/projects/1349/wiki/general_data_description_and_access/}
#' (accessed on January 15, 2016).
#' \url{http://ecocast.arc.nasa.gov/data/pub/gimms/3g.v0/00READMEgeo.txt}
#' (accessed on January 15, 2016).
#'
#' @examples
#' \dontrun{
#' ## Download sample data
#' gimms_dir <- paste0(getwd(), "/data")
#'
#' gimms_files <- downloadGimms(x = as.Date("2000-01-01"),
#'                              y = as.Date("2000-06-30"), dsn = gimms_dir)
#'
#' ## Rasterize files
#' gimms_raster <- rasterizeGimms(x = gimms_files, remove_header = TRUE)
#'
#' plot(gimms_raster)
#' }
#'
#' @export rasterizeGimms
#' @name rasterizeGimms
rasterizeGimms <-   function(x,
                             flag = FALSE,
                             header = NULL,
                             water2na = TRUE,
                             nodata2na = TRUE,
                             scaling = TRUE,
                             remove_header = FALSE,
                             cores = 1L,
                             filename = '',
                             ...) {

  ## stop if 'x' and 'filename' are of unequal length
  if (length(x) != length(filename) & all(nchar(filename)) > 0)
    stop("Parameters 'x' and 'filename' are of unequal length.")

  ## check 'cores'
  cores <- checkCores(cores)


  ### single core --------------------------------------------------------------

  if (cores == 1L | nchar(filename[1]) == 0) {

    ## loop over files in 'x'
    ls_rst <- lapply(1:length(x), function(i) {

      # import binary data as 'RasterLayer' and flip
      if (is.null(header))
        header <- createHeader(x[i])

      rst <- raster::raster(x[i])
      rst <- raster::t(rst)

      # set extent and projection
      ext <- raster::extent(c(-180, 180, -90, 90))
      rst <- raster::setExtent(rst, ext)
      raster::projection(rst) <- "+init=epsg:4326"

      # discard 'water-mask' pixels (optional)
      if (water2na)
        rst[rst[] == -10000] <- NA

      # discard 'nodata-mask' pixels (optional)
      if (nodata2na)
        rst[rst[] == -5000] <- NA

      # scale values (optional)
      if (scaling) {
        # create flags
        if (flag)
          rst <- rst - floor(rst/10) * 10 + 1
        # create ndvi
        else
          rst <- floor(rst/10) / 1000
      }

      # store output (optional)
      if (length(filename) >= i & nchar(filename[i]) > 0)
        rst <- raster::writeRaster(rst, filename = filename[i], ...)

      # if not otherwise specified, remove temporary header file
      if (remove_header)
        file.remove(header)

      # return raster
      return(rst)
    })


    ### multi-core -------------------------------------------------------------

  } else {

    ## initialize cluster
    cl <- parallel::makePSOCKcluster(cores)
    doParallel::registerDoParallel(cl)

    ## backup initial header specification
    header_init <- header

    ## loop over files in 'x'
    i <- 1
    ls_rst <- foreach::foreach(i = 1:length(x)) %dopar% {

      # import binary data as 'RasterLayer' and flip
      if (is.null(header))
        header <- createHeader(x[i])

      rst <- raster::raster(x[i])
      rst <- raster::t(rst)

      # set extent and projection
      ext <- raster::extent(c(-180, 180, -90, 90))
      rst <- raster::setExtent(rst, ext)
      raster::projection(rst) <- "+init=epsg:4326"

      # discard 'water-mask' pixels (optional)
      if (water2na)
        rst[rst[] == -10000] <- NA

      # discard 'nodata-mask' pixels (optional)
      if (nodata2na)
        rst[rst[] == -5000] <- NA

      # scale values (optional)
      if (scaling) {
        # create flags
        if (flag)
          rst <- rst - floor(rst/10) * 10 + 1
        # create ndvi
        else
          rst <- floor(rst/10) / 1000
      }

      # store output (optional)
      if (length(filename) >= i & nchar(filename[i]) > 0)
        rst <- raster::writeRaster(rst, filename = filename[i], ...)

      # if not otherwise specified, remove temporary header file
      if (remove_header)
        file.remove(header)

      # restore initial header specification
      header <- header_init

      # return raster
      return(rst)
    }

    ## deregister parallel backend
    parallel::stopCluster(cl)
  }

  ## return 'RasterLayer' (single file in 'x') or 'RasterStack'
  if (length(ls_rst) == 1) {
    return(ls_rst[[1]])
  } else {
    return(raster::stack(ls_rst))
  }
}
