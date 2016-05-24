#' Plot regression coefficients for EEM data
#' 
#' Plot regression coefficients for EEM data
#' 
#' @param x output variable from \code{\link[pls]{plsr}} function
#' @inheritParams drawEEM
#' @param ... (optional) arguments for \code{\link[EEM]{drawEEM}} and \code{\link[graphics]{filled.contour}} 
#' 
#' @return A figure is returned on the graphic device
#' 
#' @examples
#' data(gluten)
#' gluten_uf <- unfold(gluten) # unfold list into matrix
#' 
#' # delete columns with NA values
#' index <- colSums(is.na(gluten_uf)) == 0
#' gluten_uf <- gluten_uf[, index]
#' gluten_ratio <- as.numeric(names(gluten))
#' 
#' require(pls)
#' model <- plsr(gluten_ratio ~ gluten_uf, ncomp = 3) 
#' plotReg(model) 
#'  
#' @export

plotReg <- function(x, ncomp = NULL, ...){
    
    # get regression coefficient
    x <- getReg(x)
    
    # get ncomp if not provided
    if (is.null(ncomp)) ncomp = x$ncomp
    
    # plot
    drawEEM(x, ncomp, ...)
    
}
