
#' @name GetClusteringResults
#' @title Get clustering step result
#' @description \code{GetClusteringResults} returns the results of the clustering procedure \code{RunClustering}
#'
#' @seealso \code{\link{RunDenoising}}, \code{\link{RunClustering}}
#' 
#' @param data.array a (2D or 3D)+T array containing the original dynamic
#' sequence of images (the dataset). The last dimension is the time.
#' @param res.listdenois the list resulting from the \code{\link{RunDenoising}} procedure applied to \code{data.array}. This parameter may be replaced by the component \code{info.den} of the former.
#' @param res.cluster the list resulting from a call to \code{\link{RunClustering}}
#' @return a list containing two components \code{clust.array} and \code{clust.map}. \code{clust.array} is an array with same dimension as the original sequence \code{data.array} containing the clustered version. \code{clust.map} is an array with only spatial dimensions of \code{data.array} whose elements provide the cluster number at each location.
#'
#' @references Rozenholc, Y. and Reiss, M. (2012) \emph{Preserving time structures while denoising a dynamical image}, Mathematical Methods for Signal and Image Analysis and Representation (Chapter 12),  Florack, L. and Duits, R. and Jongbloed, G. and van~Lieshout, M.-C. and Davies, L. Ed., Springer-Verlag, Berlin
#' @references Lieury, T. and Pouzat, C. and Rozenholc, Y. (submitted) \emph{Spatial denoising and clustering of dynamical image sequence: application to DCE imaging in medicine and calcium imaging in neurons}  
#' @author Tiffany Lieury, Christophe Pouzat, Yves Rozenholc
#' @example ../R/Clustering-example.R    
#' @export
GetClusteringResults <- function
(data.array,         
 ## the kD+T array dataset (the last dimension is the time)
 res.listdenois,
 ## a list containing the results of the RunDenoising procedure
 res.cluster
 ## a list containing the results of the RunClustering procedure
 ){
    
    cluster.map = array(dim=dim(data.array)[-length(dim(data.array))]) 
    
    ## build the denoised kD+T array
    for (i in 1:length(res.cluster$clusters)) {
          
        for (j in res.cluster$clusters[[i]]) {
            ## retrieve information
            tmp = res.listdenois$info.den[[j]]
            
            if (is.null(tmp)) next
            
            ## get the spatial coordinates of the center
            mil = tmp$Cx
            
            ## replace the dynamics in data.array by its denoised version 
            eval(parse(text=paste('data.array[',paste(mil,collapse=','),',] <- res.cluster$centers[,',i,']',sep='')))
            
            ## add the cluster number to cluster map
            eval(parse(text=paste('cluster.map[',paste(mil,collapse=','),'] <- ',i,sep='')))
        }
    }
    
    ## return the kD+T array of the classified versions and a map of the cluster
    list(clust.array=data.array,clust.map=cluster.map,clust.center=res.cluster$centers)
}
