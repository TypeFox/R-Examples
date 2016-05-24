#' Delete scattering rays
#'
#' This function deletes three regions that are not related to fluorescence emission: 
#' (1) regions where emission wavelength is shorten than excitation light (Em <= Ex), 
#' (2) scattering rays and their second order light,
#' (3) regions above second-order scattering (EM >= 2*EX)
#'
#' @param EEM A list containing EEM data as created by \code{\link{readEEM}} function.
#' @param rep (optional) Regions to be deleted are to be replaced with \code{rep}: 0 or NA
#' @param first (optional) Width of region to be deleted for first order scattering rays [nm]
#' @param second (optional) Width of region to be deleted for second order scattering rays [nm]
#'
#' @return A list similar to input \code{EEM} is returned but with all scattering rays deleted. 
#'
#' @examples
#' data(applejuice)
#' drawEEM(delScattering2(applejuice, NA), 1)
#' 
#' @keywords scattering
#'
#' @export

delScattering2 <-
    function(EEM, rep = 0, first = 30, second = 40){

        # make sure that all EEMs has the same dimension
        dimMat <- sapply(EEM , dim)
        if (sum(!apply(dimMat, 2, function (x) identical(dimMat[,1], x))) > 0){
            stop("Dimensions do not match. Please check your data.")
        }
        
        # get Ex and Em value
        Ex <- as.numeric(colnames(EEM[[1]]))
        Em <- as.numeric(rownames(EEM[[1]]))
        numEx <- length(Ex)
        numEm <- length(Em)
        
        # expand Ex and Em to full grid
        Ex_grid <- t(matrix(rep(Ex, numEm), numEx, numEm))
        Em_grid <- matrix(rep(Em, numEx), numEm, numEx)
        
        # set index for regions below first-order scattering to be deleted (Em <= Ex)
        delIndex <- Ex_grid >= Em_grid
        
        # set index for regions that fell into 1st~2th scattering
        for (i in 1:2){
            increment <- switch(i, first, second)
            plusInd <- i * Ex_grid + increment > Em_grid
            minusInd <- i * Ex_grid - increment < Em_grid
            tempInd <- plusInd & minusInd
            delIndex <- delIndex | tempInd
        }
        
        # set index for regions above second-order scattering to be deleted (EM >= 2*Ex)
        delIndex <- delIndex | (Em_grid >= 2*Ex_grid)
        
        # number of samples and variables
        numSamp <- length(EEM)
        
        # replace targets with rep
        EEM_delS <- EEM
        for (i in 1:numSamp){
            EEM_delS[[i]][delIndex] <- rep
        }
        return(EEM_delS)
    }
