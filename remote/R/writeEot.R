if (!isGeneric('writeEot')) {
  setGeneric('writeEot', function(x, ...)
    standardGeneric('writeEot')) 
}

#' Write Eot* objects to disk
#'  
#' @description
#' Write Eot* objects to disk. This is merely a wrapper around 
#' \link{writeRaster} so see respective help section for details.
#' 
#' @param x an Eot* object
#' @param path.out the path to the folder to write the files to
#' @param prefix a prefix to be added to the file names (see Details)
#' @param overwrite see \link{writeRaster}. 
#' Defaults to \code{TRUE} in \code{writeEot} 
#' @param ... further arguments passed to \link{writeRaster}
#' 
#' @details 
#' \code{writeEot} will write the results of either an EotMode or an EotStack
#' to disk. For each mode the following files will be written:
#' 
#' \itemize{
#' \item \emph{pred_r} - the \emph{RasterLayer} of the correlation coefficients 
#' between the base point and each pixel of the predictor domain
#' \item \emph{pred_rsq} - as above but for the coefficient of determination
#' \item \emph{pred_rsq_sums} - as above but for the sums of coefficient of determination
#' \item \emph{pred_int} - the \emph{RasterLayer} of the intercept of the 
#' regression equation for each pixel of the predictor domain
#' \item \emph{pred_slp} - same as above but for the slope of the 
#' regression equation for each pixel of the predictor domain
#' \item \emph{pred_p} - the \emph{RasterLayer} of the significance (p-value) 
#' of the the regression equation for each pixel of the predictor domain
#' \item \emph{pred_resid} - the \emph{RasterBrick} of the reduced data 
#' for the predictor domain
#' }
#' 
#' Apart from \emph{pred_rsq_sums}, all these files are also created for 
#' the response domain as \emph{resp_*}. These will be pasted together
#' with the prefix & the respective mode so that the file names will 
#' look like, e.g.:
#' 
#' \emph{prefix_mode_n_pred_r.grd}
#' 
#' for the \emph{RasterLayer} of the predictor correlation coefficient 
#' of mode n using the standard \emph{raster} file type (.grd).
#' 
#' @seealso \link{writeRaster}
#' 
#' @examples
#' data(vdendool)
#' 
#' nh_modes <- eot(x = vdendool, y = NULL, n = 2, 
#'                 reduce.both = FALSE, standardised = FALSE, 
#'                 verbose = TRUE)
#' 
#' ## write the complete EotStack
#' writeEot(nh_modes, prefix = "vdendool")
#' 
#' ## write only one EotMode
#' writeEot(nh_modes[[2]], prefix = "vdendool")
#'
#' @export 
#' @name writeEot
#' @rdname writeEot 
#' @aliases writeEot,EotMode-method

# set methods -------------------------------------------------------------

setMethod('writeEot', signature(x = 'EotMode'), 
          function(x, 
                   path.out = ".",
                   prefix = "remote",
                   overwrite = TRUE,
                   ...) { 
            
            out.name <- lapply(slotNames(x)[7:19], 
                               function(i) {
                                 paste(prefix, "mode", sprintf("%02.f", 
                                                               x@mode), 
                                       i, sep = "_")
                               })
            
            out.object <- lapply(slotNames(x)[7:19], function(j) {
              slot(x, j)
            })
            
            a <- b <- NULL
            
            foreach(a = unlist(out.object), 
                    b = unlist(out.name)) %do% { 
                      raster::writeRaster(a, paste(path.out, b, sep = "/"), 
                                          overwrite = overwrite, ...)
                    }
          }
)

#' @describeIn writeEot

setMethod('writeEot', signature(x = 'EotStack'), 
          function(x, 
                   path.out = ".",
                   prefix,
                   ...) {
            
            for (i in seq(nmodes(x))) {
              writeEot(x[[i]], 
                       path.out = path.out,
                       prefix = prefix,
                       ...)
            }
          }
)