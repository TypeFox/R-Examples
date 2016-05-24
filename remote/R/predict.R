# set methods -------------------------------------------------------------
if ( !isGeneric('predict') ) {
  setGeneric('predict', function(object, ...)
    standardGeneric('predict'))
}

#' EOT based spatial prediction
#'   
#' @description
#' Make spatial predictions using the fitted model returned by \code{eot}.
#' A (user-defined) set of \emph{n} modes will be used to model the outcome 
#' using the identified link functions of the respective modes which are
#' added together to produce the final prediction.
#' 
#' @param object an Eot* object
#' @param newdata the data to be used as predictor
#' @param n the number of modes to be used for the prediction.
#' See \code{\link{nXplain}} for calculating the number of modes based 
#' on their explnatory power.
#' @param ... further arguments to be passed to \link{calc}
#' 
#' @return
#' a \emph{RasterStack} of \code{nlayers(newdata)}
#' 
#' @examples
#' ### not very useful, but highlights the workflow
#' data(pacificSST)
#' data(australiaGPCP)
#' 
#' ## train data using eot()
#' train <- eot(x = pacificSST[[1:10]], 
#'              y = australiaGPCP[[1:10]], 
#'              n = 1)
#' 
#' ## predict using identified model
#' pred <- predict(train, 
#'                 newdata = pacificSST[[11:20]], 
#'                 n = 1)
#' 
#' ## compare results
#' opar <- par(mfrow = c(1,2))
#' plot(australiaGPCP[[13]], main = "original", zlim = c(0, 10))
#' plot(pred[[3]], main = "predicted", zlim = c(0, 10))
#' par(opar)
#' 
#' @export
#' @name predict
#' @rdname predict
#' @aliases predict,EotStack-method

setMethod('predict', signature(object = 'EotStack'), 
          function(object, 
                   newdata, 
                   n = 1, 
                   ...) {
            
            ### extract identified EOT (@cell_bp) 
            ts.modes <- sapply(seq(n), function(i) {
              newdata[object[[i]]@cell_bp]
            })
            
            ### prediction using claculated intercept, slope and values
            pred.stck <- lapply(seq(raster::nlayers(newdata)), function(i) {
              raster::stack(lapply(seq(ncol(ts.modes)), function(k) {
                object[[k]]@int_response + 
                  object[[k]]@slp_response * ts.modes[i, k]
              }))
            })
            
            ### summate prediction for each mode at each time step
            pred <- raster::stack(lapply(seq(nrow(ts.modes)), function(i) {
              raster::calc(pred.stck[[i]], fun = sum, ...)
            }))
            
            return(pred)
          }
)

#' @describeIn predict

setMethod('predict', signature(object = 'EotMode'), 
          function(object, 
                   newdata, 
                   n = 1, 
                   ...) {
            
            ### extract identified EOT (@cell_bp) 
            ts.modes <- sapply(seq(n), function(i) {
              newdata[object@cell_bp]
            })
            
            ### prediction using claculated intercept, slope and values
            pred.stck <- lapply(seq(raster::nlayers(newdata)), function(i) {
              raster::stack(lapply(seq(ncol(ts.modes)), function(k) {
                object@int_response + 
                  object@slp_response * ts.modes[i, k]
              }))
            })
            
            ### summate prediction for each mode at each time step
            pred <- stack(lapply(seq(nrow(ts.modes)), function(i) {
              raster::calc(pred.stck[[i]], fun = sum, ...)
            }))
            
            return(pred)
          }
)