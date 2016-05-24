#' Cut a Minimax Linkage Tree To Get a Clustering
#' 
#' Cuts a minimax linkage tree to get one of n - 1 clusterings.  Works like
#' \code{\link{cutree}} except also returns the prototypes of the resulting
#' clustering.
#' 
#' Given a minimax linkage hierarchical clustering, this function cuts the tree
#' at a given height or so that a specified number of clusters is created.  It
#' returns both the indices of the prototypes and their locations.  This latter
#' information is useful for plotting a dendrogram with prototypes (see
#' \code{\link{plotwithprototypes}}).  As with \code{cutree}, if both k and h
#' are given, h is ignored. Unlike \code{cutree}, in current version k and h
#' cannot be vectors.
#' 
#' @param hc an object returned by \code{protoclust}
#' @param k the number of clusters desired
#' @param h the height at which to cut the tree
#' @return A list corresponding to the clustering from cutting tree:
#' \item{cl}{vector of cluster memberships} \item{protos}{vector of prototype
#' indices corresponding to the k clusters created.  \code{protos[i]} gives the
#' index of the prototype for all elements with \code{cl==i}}
#' \item{imerge}{vector describing the nodes where prototypes occur. We use the
#' naming convention of the \code{merge} matrix in \code{hclust}: if
#' \code{imerge[i]} is positive, it is the interior node (counting from the
#' bottom) of the cluster with elements \code{which(cl==i)}; if
#' \code{imerge[i]} is negative, then this is a singleton cluster with a leaf
#' as prototype.}
#' @author Jacob Bien and Rob Tibshirani
#' @seealso \code{\link{protoclust}}, \code{\link{cutree}},
#' \code{\link{plotwithprototypes}}
#' @references Bien, J., and Tibshirani, R. (2011), "Hierarchical Clustering
#' with Prototypes via Minimax Linkage," accepted for publication in \emph{The
#' Journal of the American Statistical Association}, DOI:
#' 10.1198/jasa.2011.tm10183.
#' @keywords cluster
#' @examples
#' 
#' # generate some data:
#' set.seed(1)
#' n <- 100
#' p <- 2
#' x <- matrix(rnorm(n * p), n, p)
#' rownames(x) <- paste("A", 1:n, sep="")
#' d <- dist(x)
#' 
#' # perform minimax linkage clustering:
#' hc <- protoclust(d)
#' 
#' # cut the tree to yield a 10-cluster clustering:
#' k <- 10 # number of clusters
#' cut <- protocut(hc, k=k)
#' h <- hc$height[n - k]
#' 
#' # plot dendrogram (and show cut):
#' plotwithprototypes(hc, imerge=cut$imerge, col=2)
#' abline(h=h, lty=2)
#' 
#' # get the prototype assigned to each point:
#' pr <- cut$protos[cut$cl]
#' 
#' # find point farthest from its prototype:
#' dmat <- as.matrix(d)
#' ifar <- which.max(dmat[cbind(1:n, pr[1:n])])
#' 
#' # note that this distance is exactly h:
#' stopifnot(dmat[ifar, pr[ifar]] == h)
#' 
#' # since this is a 2d example, make 2d display:
#' plot(x, type="n")
#' points(x, pch=20, col="lightblue")
#' lines(rbind(x[ifar, ], x[pr[ifar], ]), col=3)
#' points(x[cut$protos, ], pch=20, col="red")
#' text(x[cut$protos, ], labels=hc$labels[cut$protos], pch=19)
#' tt <- seq(0, 2 * pi, length=100)
#' for (i in cut$protos) {
#'   lines(x[i, 1] + h * cos(tt), x[i, 2] + h * sin(tt))
#' }
#' 
#' @export protocut
protocut <- function(hc, k=NULL, h=NULL) {
  if (!inherits(hc, "protoclust"))
    stop("Must input an object of class protoclust.")
  n <- length(hc$height) + 1
  if(is.null(k)) {
    if(is.null(h))
      stop("Must input either k or h.")
    if (length(h) > 1) {
      warning("Ignoring all but first element of h.")
      h <- h[1]
    }
    k <- n - sum(hc$height <= h)
  }
  else if(!is.null(h))
    warning("Ignoring h since k is provided.")
  if (length(k) > 1) {
    warning("Ignoring all but first element of h.")
    k <- k[1]
  }
  if(round(k) != k) {
    k <- round(k)
    cat("Rounding k to", k, ".\n")
  }
  if(k > n / 2) {
    # go bottom up
    iprot <- rep(FALSE, n - 1)
    leafprot <- rep(TRUE, n)
    if(k < n) {
      for(i in seq(n - k)) {
        for (j in 1:2) {
          if(hc$merge[i, j] < 0) {
            leafprot[-hc$merge[i, j]] <- FALSE
          }
          if(hc$merge[i, j] > 0) {
            iprot[hc$merge[i, j]] <- FALSE
          }
        }
        iprot[i] <- TRUE
      }
    }
  }
  else {
    # for k < n/2, going from top will be faster.
    iprot <- rep(FALSE, n - 1) # indices of hc$protos that we will use.
    iprot[n - 1] <- TRUE
    leafprot <- rep(FALSE, n)
    if(k > 1) {
      for(i in seq(k - 1)) {
        iprot[n - i] <- FALSE
        m <- hc$merge[n - i, ]
        for(j in 1:2) {
          if(m[j] < 0)
            leafprot[-m[j]] <- TRUE
          else
            iprot[m[j]] <- TRUE
        }
      }
    }     
  }
  leafprot <- which(leafprot)
  imerge <- c(which(iprot != 0), -leafprot)
  protos <- c(hc$protos[iprot], leafprot)
  cl <- cutree(hc, k=k)
  o <- order(cl[protos])# reorder protos in order of class
  list(cl=cl, protos=protos[o], imerge=imerge[o])
}
