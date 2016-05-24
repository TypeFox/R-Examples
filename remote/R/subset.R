if (!isGeneric('subset')) {
  setGeneric('subset', function(x, ...)
    standardGeneric('subset')) 
}

#' Subset modes in EotStacks
#' 
#' @description
#' Extract a set of modes from an EotStack
#' 
#' @param x EotStack to be subset
#' @param subset integer or character. The modes to ectract (either by
#' integer or by their names)
#' @param drop if \code{TRUE} a single mode will be returned as an EotMode
#' @param ... currently not used
#' 
#' @return
#' an Eot* object
#' 
#' @examples
#' data(vdendool)
#' 
#' nh_modes <- eot(x = vdendool, y = NULL, n = 3, 
#'                 reduce.both = FALSE, standardised = FALSE, 
#'                 verbose = TRUE)
#'                 
#' subs <- subset(nh_modes, 2:3) # is the same as
#' subs <- nh_modes[[2:3]]
#' 
#' ## effect of 'drop=FALSE' when selecting a single layer
#' subs <- subset(nh_modes, 2)
#' class(subs)
#' subs <- subset(nh_modes, 2, drop = TRUE)
#' class(subs)
#' 
#' @export
#' @name subset
#' @rdname subset
#' @aliases subset,EotStack-method
#' @aliases subset,ANY-method

setMethod('subset', signature(x = 'EotStack'), 
          function(x, subset, drop = FALSE, ...) {
            if (is.character(subset)) {
              i <- na.omit(match(subset, names(x)))
              if (length(i) == 0) {
                stop('invalid mode names')
              } else if (length(i) < length(subset)) {
                warning('invalid mode names omitted')
              }
              subset <- i
            }
            subset <- as.integer(subset)
            if (! all(subset %in% 1:nmodes(x))) {
              stop('not a valid subset')
            }
            if (length(subset) == 1 & drop) {
              x <- x@modes[[subset]]
            } else {
              x@modes <- x@modes[subset]
            }
            return(x)
          }
)

#' @describeIn subset
#' @param i number of EotMode to be subset

setMethod("[[", signature(x = "EotStack"), 
          function(x, i) {
            subset(x, i, drop = TRUE)
          }
)
