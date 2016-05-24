#' @export thin.algorithm
#' @title Implements random spatial thinning algorithm
#' 
#' @description \code{thin.algorithm} implements a randomization approach to
#' spatially thinning species occurence data. This function is the algorithm underlying
#' the \code{\link{thin}} function.
#' 
#' @param rec.df.orig A data frame of long/lat points for each presence record. The 
#' data.frame should be a two-column data frame, one column of long and one of 
#' lat
#' @param thin.par Thinning parameter - the distance (in kilometers) that you want
#' records to be separated by.
#' @param reps The number of times to repete the thinning process. Given the random
#' process of removing nearest-neighbors there should be 'rep' number of different
#' sets of coordinates.
#' @return reduced.rec.dfs: A list object of length 'rep'. Each list element is a different
#' data.frame of spatially thinned presence records.
#' 


thin.algorithm <- function( rec.df.orig, thin.par, reps ) {

  ## Create empty list object to store thinned occurrence datasets
  reduced.rec.dfs <- list()
  
  for ( Rep in 1:reps ){
    
    rec.df <- rec.df.orig
    
    ## Calculate distance matrix using fields::rdist.earth
    ## ***
    ## This function returns distances, correcting for change in distance
    ## of degree of longitude based on distance from equator
    DistMat <- rdist.earth( x1=rec.df, miles=FALSE )
    
    ## Set diagonal elements of distance matrix to NAs, as these are all 
    ## 0 in this case.
    diag(DistMat) <- NA
    
    ## Perform while loop based on two criteria
    ## 1. The minimum distance between to occurences is less than the 
    ##    thinning parameter 
    ## 2. The number of rows in the data.frame is greater than 1
    while( min( DistMat, na.rm=TRUE ) < thin.par & nrow( rec.df ) > 1 ) {
      
      ## Find array indices for all occurrences that are less than thin.par
      ## away from another occurence, and return the distance matrix **row** 
      ## values for these occurrences
      CloseRecs <- which( DistMat < thin.par, arr.ind=TRUE )[ , 1]
      
      ## Multiple steps in this one line:
      ## ***
      ## a. For each row (occurrence) in the distance matrix, determine how
      ##    many other occurrences are within thin.par away.
      ## b. Determine which occurence(s) have the greatest number of occurrences
      ##    within thin.par distance.
      ## c. Identify these occurences by name and convert them into numeric values
      ##    (to be used as indexs of occurrences to remove).
      RemoveRec <- as.numeric( names( which( table( CloseRecs ) == max( table( CloseRecs ) ) ) ) )
      
      ## If there are more than one occurrences meeting the above
      ## conditions, choose one to remove at random.
      if( length( RemoveRec ) > 1 ) {
        RemoveRec <- sample( RemoveRec, 1 )
      }
      
      ## Remove a single occurrence that has the most occurrences within 
      ## thin.par distance from the dataset
      rec.df <-rec.df[ -RemoveRec, ]
      
      ## Remove the occurence from the distance matrix
      DistMat <- DistMat[ -RemoveRec, -RemoveRec ]
      
      ## Break out of while loop if there is only one record left
      if( length( DistMat ) == 1 ){ break }
      
    }
    
    ## Save the thinned dataset to reduced.rec.dfs data.frame
    colnames(rec.df) <- c("Longitude", "Latitude")
    reduced.rec.dfs[[Rep]] <- rec.df
  }
  return( reduced.rec.dfs )
}


## ******************************************************************** ##
## Below is legacy code that is temporarly being stored in this file.
## After a round of peer-review, we recieved a suggestion that in how
## the algorithm works that fundamentally changed the alogirthm and
## increased its efficiency by orders of magnitude.
## The code below remains so that MA-L can incorporate comments
## appropriately to the new code.
## ******************************************************************** ##
#   ##
#   ## General Algorithm:
#   ## ******************
#   ## While at least one member of the nearest neighbor vector is less than 10 km
#   ## remove one of the points associated with nearest neighbor less than 10km randomly
#   ##
#   ## A few more details:
#   ## 1. Define thinning parameter
#   ## 2. Find all points with nearest neighbor closer than thinning parameter
#   ## 3. Randomly delete one of these points
#   ## 4. Continue until no points have nearest neighbor closer than thinning parameter
#   ## 5. Repeat many times 
#   ##      
#   
#   reduced.rec.dfs <- list()
#   
#   # Value for which neighber is too close - i.e. thinning metric
#   print( paste( "Thinning data using thin.par: ", thin.par ) )
#   
#   for ( Rep in 1:reps ){
#     #print( paste( "Repetition: ", Rep ) ) ### Debug Line - uncomment to print repetition number
#     rec.df <- rec.df.orig
# 
#     # Make distance matrix
#     DistMat <- rdist.earth(x1=rec.df,miles=FALSE)
#     # print(DistMat) ### DEBUG LINE - warning this could be very large!
#     # Replace diagonal (which is all zeroes) with NAs
#     diag(DistMat) <- NA
#     # For each record, calculate nearest neighbor
#     NearNeighbor <- apply( DistMat, MARGIN=2, FUN=min, na.rm=TRUE )
#     MinNearNeighbor <- min(NearNeighbor)
#     # print(MinNearNeighbor) ### DEBUG LINE
#     
#     while( MinNearNeighbor<thin.par & nrow(rec.df)>1 ){
#       
#       # Determine records with a minimum nearest neighbor distance less than thin.par 
#       CloseRecs <- which( NearNeighbor < thin.par )
#       # print( paste( 'Length of CloseRecs: ',length(CloseRecs)) ) ### DEBUG LINE
#       
#       # Choose one record to remove randomly
#       RemoveRec <- sample(CloseRecs,1)
#       # print( paste( 'RemoveRec: ',RemoveRec) ) ### DEBUG LINE
#       
#       # Remove random record from the complete records data.frame rec.df
#       rec.df <-rec.df[ -RemoveRec, ]
#       # print( paste( 'Number of remaining records: ', nrow(rec.df) )) ### DEBUG LINE
#       
#       # Remove random record from the distance matrix DistMat
#       DistMat <- DistMat[ -RemoveRec, -RemoveRec ]
# 
#       # Calculate NEW nearest neighbors **UNLESS* there is only one location
#       # remaining. 
#       if ( nrow(rec.df)>1 ){
#         NearNeighbor <- apply( DistMat, MARGIN=2, FUN=min, na.rm=TRUE )
#         MinNearNeighbor <- min(NearNeighbor)
#         # print(paste( 'MinNearNeighbor inside loop: ', MinNearNeighbor) ) DEBUG LINE
#       }
#     }
#     reduced.rec.dfs[[Rep]] <- rec.df
#   }
#   return( reduced.rec.dfs )
# }