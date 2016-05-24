#' Nearest Neighbor Distances for Morphological Disparity Studies
#'
#' This is a simple function for obtaining nearest neighbor distance
#' from a symmetric pair-wises distance matrix, assumed here to be 
#' dissimilarities between pairs of taxa. Per-species NND is returned
#' rather than a mean or other summary value.

#' @details
#' This function is mainly included here for pedagogical (teaching) purposes.
#' NND is so simple to calculate, users are urged to write their own functions
#' for primary research purposes.
#'
#' Typically, the \emph{mean} NND for a group is reported and used to compare different
#' groupings of taxa (such as different time intervals, or different clades). Bootstrapping
#' should be used to generate confidence intervals.

#' @param distMat A symmetric, square pair-wise distance matrix, assumed to be a dissimilarity
#' matrix with a zero diagonal. Can be a 'dist' object rather than a numerical matrix.
#' Taxon labels can be applied to the rows and columns (or as labels if a 'dist' object)
#' and will be used to name the resulting output.

#' @return
#' Returns a vector of the nearest neighbor distance for each unit (taxon) in the 
#' pair-wise distance matrix, named with the labels from the input distance matrix.

#' @aliases nearestNeighborDist

#' @seealso
#' For the example dataset used in examples, see \code{\link{graptDisparity}}

#' @author David W. Bapst

#' @references
#' Bapst, D. W., P. C. Bullock, M. J. Melchin, H. D. Sheets, 
#' and C. E. Mitchell. 2012. Graptoloid diversity and disparity 
#' became decoupled during the Ordovician mass extinction. \emph{Proceedings 
#' of the National Academy of Sciences} 109(9):3428-3433.
#'
#' Ciampaglio, C. N., M. Kemp, and D. W. McShea. 2001. Detecting 
#' changes in morphospace occupation patterns in the fossil record: 
#' characterization and analysis of measures of disparity. \emph{Paleobiology}
#' 27(4):695-715.
#' 
#' Foote, M. 1990. Nearest-neighbor analysis of trilobite morphospace.
#' \emph{Systematic Zoology} 39:371-382.


#' @examples
#' #example using graptolite disparity data from Bapst et al. 2012
#' 
#' #load data
#' data(graptDisparity)
#' 
#' #calculate mean NND
#' NND<-nearestNeighborDist(graptDistMat)
#' mean(NND)
#' 
#' #calculate NND for different groups
#' 
#' #group (clade/paraclade) coding
#' groupID <- graptCharMatrix[,54]+1
#' 
#' groupNND<-numeric(7)
#' names(groupNND)<-c("Normalo.","Monogr.","Climaco.",
#'    "Dicrano.","Lasiogr.","Diplogr.","Retiol.")
#' for(i in unique(groupID)){
#'    groupNND[i]<-mean(nearestNeighborDist(
#'       graptDistMat[groupID==i,groupID==i]))
#'    }
#' groupNND
#' 
#' #the paraphyletic Normalograptids that survived the HME are most clustered
#'    #but this looks at all the species at once
#'    #and doesn't look for the nearest *co-extant* neighbor!
#'    #need to bring in temporal info to test that
#' 

#' @name nearestNeighborDist
#' @rdname nearestNeighborDist
#' @export nearestNeighborDist
nearestNeighborDist<-function(distMat){
   #returns a vector of NNDs for each taxon in a distance matrix
   #this function is included in paleotree mainly for pedagogical use
   if(!inherits(distMat,"dist")){
      if(is.matrix(distMat)){
         if(all(diag(distMat)!=0)){
            stop("Diagonal is nonzero, may be a similarity matrix, not a distance matrix")}
         if(!isSymmetric(distMat)){stop("Not a Symmetric Distance Matrix?")}
         distM<-distMat
      }else{
         stop("Not a matrix, not a 'dist' object, what is it?")
      }
   }else{
      distM<-is.matrix(distMat)
      }
   NND<-sapply(1:nrow(distMat),function(x) min(distMat[x,-x]))
   names(NND)<-labels(distM)[[1]]
   return(NND) #return as a per-taxon 
   }
   
