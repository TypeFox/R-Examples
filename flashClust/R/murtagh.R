# Code by F. Murthagh, http://astro.u-strasbg.fr/~fmurtagh/mda-sw/splus
# modified by Peter Langfelder to make it compaticle with R's standard hclust

flashClust <- function(d, method="complete",  members = NULL)
{
  hclust(d, method, members)
}

hclust <- function(d, method="complete",  members = NULL)
{

# Hierarchical clustering, on raw input data; we will use Euclidean distance.
# A range of criteria are supported; also there is a storage-economic option.
# Author: F. Murtagh, May 1992

 METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
        "median", "centroid")
 method <- pmatch(method, METHODS)

 if (is.na(method)) 
        stop("Invalid clustering method")
    if (method == -1) 
        stop("Ambiguous clustering method")
 n = attr(d, "Size")
 len = length(d);
 if (len!=(n*(n-1)/2))
   stop("Distance structure appears invalid.");

 if (n==1 || len==0)
   stop("The distance structure is empty.");

 if (is.null(members)) 
 {
        members <- rep(1, n)
 }
 else if (length(members) != n) 
        stop("invalid length of members")

# We choose the general routine, `hc', which
# caters for 7 criteria, using a half dissimilarity matrix; (BTW, this uses the
# very efficient nearest neighbor chain algorithm, which makes this algorithm
# of O(n^2) computational time, and differentiates it from the less efficient
# -- i.e. O(n^3) -- implementations in all commercial statistical packages
# -- as far as I am aware -- except Clustan.)  

 hcl <- .Fortran("hc",
          n = as.integer(n),
          len = as.integer(len),
          method = as.integer(method),
          ia = integer(n),
          ib = integer(n),
          crit = double(n),
          membr = as.double(members),
          nn = integer(n),
          disnn = double(n),
          flag = logical(n),
          diss = as.double(d), 
          PACKAGE = "flashClust")

# 2nd step: interpret the information that we now have, -- seq. of aggloms., --
# as merge, height, and order lists.

#PL: not clear what this iclass is supposed to be for.

 #iclass <- matrix(0.0, n, n)
 #storage.mode(iclass) <- "integer"

 hcass <- .Fortran("hcass2",
          n = as.integer(n),
          ia = as.integer(hcl$ia),
          ib = as.integer(hcl$ib),
          order = integer(n),
          iia = integer(n),
          iib = integer(n), 
          PACKAGE = "flashClust")

 merge <- cbind(hcass$iia[1:n-1],hcass$iib[1:n-1])


 hhh <- list(merge = merge, height = hcl$crit[1:n-1], order = hcass$order,
             labels = attr(d, "Labels"), method = METHODS[method], 
             call = match.call(), dist.method = attr(d, "method"))
 class(hhh) = "hclust"


 hhh

}




