#' Extract values from other models
#' 
#' Extract values from other models
#' 
#' @param x output variable from \code{\link[stats]{prcomp}} or \code{\link[pls]{plsr}} functions
#' 
#' @return A `EEMweight` list containing title and value attributes. 
#' 
#' @examples
#' data(applejuice)
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' result <- prcomp(applejuice_uf) 
#' loading <- getLoading(result)
#' str(loading)
#' @name extract
NULL

#' @rdname extract
#' @export
getLoading <- function(x) {            
    # extract information
    xClass <- class(x)
    if (xClass == "prcomp") {
        loading <- x$rotation
        loadingNames <- rownames(loading)
        title = "Loading"
        ncomp = NULL
    } else if (xClass == "mvr") {
        loading <- x$loadings
        loadingNames <- rownames(loading)
        title = "LV"
        ncomp = x$ncomp
    } else {
        stop('Input class not supported.')
    }
    rownames(loading) <- loadingNames # EX...EM... format
    y <- list(title = title, value = loading, ncomp = ncomp) 
    class(y) <- "EEMweight"
    return(y)
}

#' @rdname extract
#' @export
getReg <- function(x) {
    
    # extract information
    if (class(x) == "mvr") {
        tmpreg <- x$coefficients
        reg <- as.matrix(tmpreg[,1,])
        regNames <- rownames(x$coefficients)
        title = "Regression coefficient"
        ncomp = x$ncomp
    } else {
        stop('Input class not supported.')
    }
    rownames(reg) <- regNames # EX...EM... format
    y <- list(title = title, value = reg, ncomp = ncomp)
    class(y) <- "EEMweight"
    return(y)
}
