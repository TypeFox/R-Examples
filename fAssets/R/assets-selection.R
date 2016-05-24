
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                   DESCRIPTION:
#  assetsSelect                Selects similar or dissimilar assets 
#  .hclustSelect               Selects due to hierarchical clustering 
#  .kmeansSelect               Selects due to k-means clustering    
################################################################################


assetsSelect <- 
    function(x, method = c("hclust", "kmeans"), control = NULL, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description: 
    #   Clusters a set of assets
    
    # Arguments:
    #   method - which algorithm should be used?
    #       hclust - Hierarchical clustering on a set of dissimilarities 
    #       kmeans - k-means clustering on a data matrix 
        
    # FUNCTION:

    # Selection:
    # do not method = match.arg(method) to allow for user specified clustering
    method <- method[1]   
     
    # Transform to matrix:
    if (class(x) == "timeSeries") {
        x <- as.matrix(x)
    }
    
    # Compose Function:
    fun <- paste(".", method, "Select", sep = "")
    FUN <- get(fun)

    # Cluster:
    ans <- FUN(x, control, ...)
    
    # Return Value:
    ans
}


################################################################################


.hclustSelect <-  
    function(x, control = NULL, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description: 
    #   Hierarchical Clustering
    
    # FUNCTION:
    
    # Method:
    if (is.null(control)) 
        control = c(measure = "euclidean", method = "complete")
    measure = control[1]
    method = control[2]
    
    # hclust:
    ans = hclust(dist(t(x), method = measure), method = method, ...)
    class(ans) = c("list", "hclust")
    
    # Return Value:
    ans
}


################################################################################


.kmeansSelect <-  
    function(x, control = NULL, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description: 
    #   kmeans Clustering
    
    # Note:
    #   centers must be specified by the user!
    
    # FUNCTION:
    
    # Method:
    if (is.null(control)) 
        control = c(centers = 5, algorithm = "Hartigan-Wong")
    centers = as.integer(control[1])
    algorithm = control[2]
    
    # kmeans:
    ans = kmeans(x = t(x), centers = centers, algorithm = algorithm, ...)
    class(ans) = c("list", "kmeans")
    
    # Return Value:
    ans
}


################################################################################