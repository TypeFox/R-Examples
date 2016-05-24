#' Plot loadings for EEM data
#' 
#' Plot loadings for EEM data
#' 
#' @param x output variable from \code{\link[stats]{prcomp}} or \code{\link[pls]{plsr}} functions
#' @inheritParams drawEEM
#' @param ... (optional) arguments for \code{\link[EEM]{drawEEM}} and \code{\link[graphics]{filled.contour}} 
#' 
#' @return A figure is returned on the graphic device
#' 
#' @examples
#' data(applejuice)
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' result <- prcomp(applejuice_uf) 
#' plotLoading(result, ncomp = 1) # plot loading of the first PC
#'  
#' @export
#' 

plotLoading <- function(x, ncomp = NULL, ...){
        
        # get loading
        x <- getLoading(x)
        
        # get ncomp if not provided
        if (is.null(ncomp)) {
            # call ncomp from EEMweight. If not present abort the function.
            ncomp = x$ncomp
            if (is.null(ncomp)) stop("ncomp must be provided")
            }

        # plot
        drawEEM(x, ncomp, ...)

    }
