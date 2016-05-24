
#' @name GetDenoisingResults
#' @title Get denoising step result
#' @description \code{GetDenoisingResults} returns the denoised version of a dynamical image sequence
#' as an array having the same dimensions as the original sequence.
#'
#' @seealso \code{\link{RunDenoising}}
#' 
#' @param data.array a (2D or 3D)+T array containing the original dynamic
#' sequence of images (the dataset). The last dimension is the time.
#' @param res.listdenois the list resulting from the \code{\link{RunDenoising}} procedure applied to \code{data.array}. This parameter may be replaced by the component \code{info.den} of the former.
#' @return an array with same dimension as \code{data.array} containing the denoised version.
#'
#' @references Rozenholc, Y. and Reiss, M. (2012) \emph{Preserving time structures while denoising a dynamical image}, Mathematical Methods for Signal and Image Analysis and Representation (Chapter 12),  Florack, L. and Duits, R. and Jongbloed, G. and van~Lieshout, M.-C. and Davies, L. Ed., Springer-Verlag, Berlin
#' @references Lieury, T. and Pouzat, C. and Rozenholc, Y. (submitted) \emph{Spatial denoising and clustering of dynamical image sequence: application to DCE imaging in medicine and calcium imaging in neurons}  
#' @author Tiffany Lieury, Christophe Pouzat, Yves Rozenholc
#' @example ../R/Denoising-example.R
#' @export
GetDenoisingResults <- function
(data.array,         
 ## the kD+T array dataset (the last dimension is the time)
 res.listdenois
 ## a list containing the results of the RunDenoising procedure
 ){
    
    ## build the denoised kD+T array
    for (i in 1:length(res.listdenois$info.den)) {
        
        ## retrieve information
        tmp = res.listdenois$info.den[[i]]
        
        if (is.null(tmp)) next
        
        ## get the spatial coordinates of the center
        mil = tmp$Cx
        
        ## replace the dynamics in data.array by its denoised version 
        eval(parse(text=paste('data.array[',paste(mil,collapse=','),',] <- tmp$Ix',sep='')))
    }
    
    ## return the kD+T array of the denoised versions
    data.array
}
