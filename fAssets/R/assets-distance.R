
# This library is free software, you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation, either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library, if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  assetsDist            Computes the distances between assets
# FUNCTION:             DESCRIPTION:
#  corDist               Returns correlation distance measure
#  kendallDist           Returns kendalls correlation distance measure        
#  spearmanDist          Returns spearmans correlation distance measure
#  mutinfoDist           Returns mutual information distance measure
# FUNCTION:             DESCRIPTION:
#  euclideanDist         Returns Euclidean distance measure
#  maximumDist           Returns maximum distance measure
#  manhattanDist         Returns Manhattan distance measure
#  canberraDist          Returns Canberra distance measure
#  binaryDist            Returns binary distance measure
#  minkowskiDist         Returns Minkowsky distance measure
# FUNCTION:             DESCRIPTION:
#  braycurtisDist        Returns Bray Curtis distance measure
#  mahalanobisDist       Returns Mahalanobis distance measure
#  jaccardDist           Returns Jaccard distance mesaure
#  sorensenDist          Returns Sorensen distance measure
################################################################################


assetsDist <-
    function(x, method="cor", ...)
{
     fun <- match.fun(paste(method, "Dist", sep=""))
     if (method == "mutinfo") {
         ans <- fun(x, ...)
     } else {
         ans <- fun(x)
     }
     
     # Return Value:
     ans  
}
     

################################################################################


corDist <-
    function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns correlation distance
    
    # Argument:
    #   x - a bivariate- or multivariate 'timeSeries' object
    
    # Example:
    #   corDist(matrix(rnorm(100), ncol=5))
    
    # FUNCTION:
    
    # Distance:
    dist <- as.dist(1-cor(x))
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


kendallDist <-
    function(x)
{ 
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Kendal's correlation distance
    
    # FUNCTION:
    
    # Distance:
    dist <- as.dist(1-cor(x, method = "kendall"))
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


spearmanDist <-
    function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Spearman's correlation distance
    
    # FUNCTION:
    
    # Distance:
    dist <- as.dist(1-cor(x, method = "spearman"))
    
    # Return Value:
    dist
}


###############################################################################
# Selected distance measures from contributed R package bioDist:


mutinfoDist <-
    function(x, nbin=10)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns mutual information distance measure
    
    # Note:
    #   borrowed from R package bioDist and slightly modified
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    nc <- ncol(x)
    nr <- nrow(x)
    clist <- vector("list", length=nr)
    for(i in 1:nr) clist[[i]] <- cut(x[i,], breaks=nbin)
    
    ppfun <- function(pp) 
    {
        pp <- pp[pp > 0]
        -sum(pp*log(pp ))
    }
    appfun <- function(x,y) 
    {
        ppfun(table(x)/nc) + ppfun(table(y)/nc) - ppfun(c(table(x, y)/nc))
    }
        
    mat <- matrix(rep(NA, nr*nr), ncol = nr)
    for(i in 1:(nr-1)) 
    {
        for(j in (i+1):nr) 
        {
            mat[i,j] <- mat[j,i]<- appfun(clist[[i]], clist[[j]])
        }
    }
    
    # Distance:
    mat <- 1 - sqrt(1 - exp(-2*mat))
    colnames(mat) <- rownames(mat) <- rownames(x)
    dist <- as.dist(mat)
        
    # Return Value:
    dist
}
    
    
################################################################################
# Selected distance functions from base R Package:
# "euclidean", "maximum", "manhattan", 
# "canberra", "binary", "minkowski" 
  
    
euclideanDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Euclidean distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- dist(x, "euclidean")
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


maximumDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns maximum distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- dist(x, "maximum")
    
    # Return Value:
    dist  
}


# ------------------------------------------------------------------------------


manhattanDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Manhattan distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- dist(x, "manhattan")
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


canberraDist <-
function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Canberra distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- dist(x, "canberra")
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


binaryDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns binary distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- dist(x, "binary")
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


minkowskiDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Minkowski distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- dist(x, "minkowski")
    
    # Return Value:
    dist
}


################################################################################
# Selected distance from contributed R package ecodist:
# "bray-curtis", "mahalanobis", "jaccard", "sorensen" 
# See builtin script: dist-distEcodist.R
    

braycurtisDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # FUNCTION:
    
    # Distance:
        x <- t(as.matrix(x))
    dist <- ecodist::bcdist(x)
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


mahalanobisDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Mahalanobis distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- ecodist::distance(x, "mahalanobis")
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


jaccardDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns Jaccard's distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- ecodist::distance(x, "jaccard")
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


sorensenDist <-
    function(x)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns difference distance measure
    
    # FUNCTION:
    
    # Distance:
    x <- t(as.matrix(x))
    dist <- ecodist::distance(x, method="sorensen")
    
    # Return Value:
    dist 
}


################################################################################

