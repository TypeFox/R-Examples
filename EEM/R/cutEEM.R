#' Cut portions of EEM
#'
#' Cut portions of EEM
#'
#' @inheritParams drawEEM
#' @param cutEX Numeric or sequential data specifying regions to be cut for excitation wavelength. 
#' Examples, 200 or 200:500 
#' @param cutEM Numeric or sequential data specifying regions to be cut for emission wavelength. 
#' Examples, 200 or 200:500 
#'
#' @return A list similar to input \code{EEM} is returned but with specified portions cut. 
#'
#' @examples
#' data(applejuice)
#' applejuice_cut <- cutEEM(applejuice, cutEX = 300:450)
#' drawEEM(applejuice_cut, 1)
#'
#' @export

cutEEM <- function(x, cutEX = NULL, cutEM = NULL) UseMethod("cutEEM", x)

#' @rdname cutEEM
#' @export
cutEEM.EEM <- function(x, cutEX = NULL, cutEM = NULL){
    
    # check that all x has the same dimension
    dimMat <- sapply(x, dim)
    if (sum(!apply(dimMat, 2, function (x) identical(dimMat[,1], x))) > 0){
        stop("Dimension do not match. Please check your data.")
    }    
    
    # prepare data 
    N <- length(x)
    rowname <- as.numeric(rownames(x[[1]])) # EM
    colname <- as.numeric(colnames(x[[1]])) # EX
    
    # check that it cannot cut through the middle
    # for EX
    if (!is.null(cutEX)){
        if (suppressWarnings(min(colname) < min(cutEX)&
                                 max(colname) > max(cutEX))){
            stop("Cannot cut through the middle.")
        }
    } 
    # for EM
    if (!is.null(cutEM)){    
        if (suppressWarnings(min(rowname) < min(cutEM)&
                                 max(rowname) > max(cutEM))) {
        stop("Cannot cut through the middle.")
        }
    }
    
    # function to find logical index to cut
    cut_index <- function(WL, cut_value){
        if (is.null(cut_value)){
            rep(FALSE, length(WL))
        } else {
            min(cut_value) <= WL & WL <= max(cut_value)
        }
    }
    
    # prepare logical value for rows remain after cutting (EM)
    rowIdx <- !cut_index(rowname, cutEM)
    
    # prepare columns remain after cutting (EX)
    colIdx <-  !cut_index(colname, cutEX)
    colSelected <- which(colIdx)
    
    # cut
    x_cut <- list()
    for (i in 1:N){
        x_cut[[i]] <- subset(x[[i]],
                             subset = rowIdx, # subset refer to rows
                             select = colIdx,  # select refer to columns
                             drop = TRUE)
    }
    names(x_cut) <- names(x)
    class(x_cut) <- "EEM"
    return(x_cut)
}

#' @rdname cutEEM
#' @export
cutEEM.EEMweight <- function(x, cutEX = NULL, cutEM = NULL){
    
    # get value 
    title <- x$title
    EEM <- x$value

    # prepare data
    varnames <- rownames(x)
    EX <- getEX(varnames)
    EM <- getEM(varnames)
    
    # check that it cannot cut through the middle
    # for EX
    if (!is.null(cutEX)){
        if (suppressWarnings(min(EX) < min(cutEX)&
                                 max(EX) > max(cutEX))){
            stop("Cannot cut through the middle.")
        }
    } 
    # for EM
    if (!is.null(cutEM)){    
        if (suppressWarnings(min(EM) < min(cutEM)&
                                 max(EM) > max(cutEM))) {
            stop("Cannot cut through the middle.")
        }
    }
    
    # function to find logical index to cut
    cut_index <- function(WL, cut_value){
        if (is.null(cut_value)){
            rep(FALSE, length(WL))
        } else {
            min(cut_value) <= WL & WL <= max(cut_value)
        }
    }
    
    # prepare logical index for columns remaining after cutting 
    EXIdx <-  !cut_index(EX, cutEX)
    EMIdx <- !cut_index(EM, cutEM)
    index <- EXIdx & EMIdx
    
    # cut
    x_cut <- x[index,]
    
    # return output
    output <- list()
    output$title <- title
    output$value <- x_cut
    class(output) <- "EEMweight"
    
    return(output)
}
