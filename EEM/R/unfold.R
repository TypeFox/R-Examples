#' Unfold EEM list into a matrix 
#' 
#' Unfold EEM list into a matrix with columns as variables (wavelength conditions) and rows as samples.
#' 
#' @param EEM a list containing EEM data as created by \code{readEEM} function.
#' @param replaceNA logical value whether to replace NA with 0
#' 
#' @return Unfolded EEM matrix where columns are wavelength condition and rows are samples
#' 
#' @examples
#' data(applejuice)
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' dim(applejuice_uf) # dimension of unfolded matrix
#' 
#' @export

unfold <- function(EEM, replaceNA = TRUE){
        
        ## check that all EEM has the same dimension
        dimMat <- sapply(EEM, dim)
        if (sum(!apply(dimMat, 2, function (x) identical(dimMat[,1], x))) > 0){
            stop("Dimension do not match. Please check your data.")
        }
        
        ## check that EX and EM wavelength are identical for all samples
        dimension_names <- lapply(EEM, dimnames)
        if (sum(!sapply(dimension_names, function (x) identical(dimension_names[[1]], x))) > 0){
            stop("Dimension names do not match. Please check your data.")
        }
        
        # get sName
        sName <- names(EEM)
        
        ## begin unfolding EEM data into a matrix with no. of samples x variables  
        EEM_uf <- sapply(1:length(EEM), function(i){
            as.numeric(EEM[[i]])
        })
        
        EEM_uf <- t(EEM_uf)
        
        # get var
        var.full <- expand.grid(em = as.numeric(rownames(EEM[[1]])), ex = as.numeric(colnames(EEM[[1]])))
        var <- paste0("EX", var.full$ex, "EM", var.full$em)
        
        colnames(EEM_uf) <- var
        rownames(EEM_uf) <- sName
        
        # replace NA with 0
        if (replaceNA) EEM_uf[is.na(EEM_uf)] <- 0
        
        return(EEM_uf)
    }
