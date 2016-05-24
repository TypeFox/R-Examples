#' Hierarchical Clustering with Prototypes: Minimax Linkage.
#' 
#' Performs minimax linkage hierarchical clustering given a set of
#' dissimilarities.  Returns an object that looks just like the output of
#' \code{hclust} except that it has an additional element containing prototype
#' indices.
#' 
#' This function provides an efficient implementation of minimax linkage
#' hierarchical clustering.  Consider two clusters G and H and their union U.
#' The minimax linkage between G and H is defined to be the radius of the
#' smallest ball that encloses all of U and that is centered at one of the
#' points in U.  If G and H are merged together, the prototype for the newly
#' formed cluster U is that enclosing ball's center.  By construction, the
#' prototype for a cluster will always be one of the objects being clustered.
#' For more on minimax linkage and how one can use prototypes to help interpret
#' a dendrogram, see
#' 
#' Bien, J., and Tibshirani, R. (2011), "Hierarchical Clustering with
#' Prototypes via Minimax Linkage," accepted for publication in \emph{The
#' Journal of the American Statistical Association}, DOI:
#' 10.1198/jasa.2011.tm10183.
#' 
#' This function has been designed to work like \code{hclust} in terms of
#' inputs and outputs; however, unlike \code{hclust}, it outputs an additional
#' element, namely a vector of length n - 1 containing the indices of
#' prototypes.  It follows \code{hclust}'s convention for making the arbitrary
#' choice of whether to put a subtree on the left or right side.
#' 
#' For cutting a minimax linkage hierarchical clustering, use
#' \code{\link{protocut}}, which works like \code{\link{cutree}} except that it
#' returns the set of prototypes in addition to the cluster assignments.
#' 
#' This function calls a C implementation of the algorithm detailed in Bien and
#' Tibshirani (2011) that is based on an algorithm described in Murtagh (1983).
#' 
#' @param d dissimilarities object.  Can be of class \code{dist} or
#' \code{matrix}
#' @param verb see verbose output?
#' @return An object of class \code{protoclust}, which is just like
#' \code{hclust} but has an additional element: \item{merge, height,
#' order}{identical to the values returned by \code{\link{hclust}}}
#' \item{protos}{a vector of length n - 1.  The i-th element is the index of
#' the prototype corresponding to the cluster formed on the i-th merge.}
#' @author Jacob Bien and Rob Tibshirani
#' @seealso \code{\link{protocut}}, \code{\link{plotwithprototypes}},
#' \code{\link{hclust}}
#' @references Bien, J., and Tibshirani, R. (2011), "Hierarchical Clustering
#' with Prototypes via Minimax Linkage," accepted for publication in \emph{The
#' Journal of the American Statistical Association}, DOI:
#' 10.1198/jasa.2011.tm10183.
#' 
#' Murtagh, F. (1983), "A Survey of Recent Advances in Hierarchical Clustering
#' Algorithms," \emph{The Computer Journal}, \bold{26}, 354--359.
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
#' @export protoclust
protoclust <- function(d, verb=FALSE) {
  if (is.matrix(d)) {
    if (nrow(d) == ncol(d)) {
      cat("converting to dist (note: ignores above diagonal)", fill=TRUE)
      d <- as.dist(d)
    }
  }
  if (class(d) != "dist")
    stop("d must be of class \"dist\" or be a square matrix.")
  nn <- length(d)
  n <- (1 + sqrt(1 + 8 * nn)) / 2
  
  out <- .C("hier",
            as.double(d),
            as.integer(n),
            as.integer(verb),
            as.integer(matrix(0, n - 1, 2)),
            as.double(rep(0, n - 1)),
            as.integer(rep(0, n)),
            as.integer(rep(0, n - 1)), PACKAGE="protoclust")[4:7]
    
  tree <- list(merge=matrix(out[[1]], n - 1, 2, byrow=TRUE),
               height=out[[2]],
               order=out[[3]],
               protos=out[[4]],
               labels=attr(d, "Labels"),
               method="minimax",
               call=match.call(),
               dist.method=attr(d, "method"))
  class(tree) <- c("protoclust", "hclust") # inherits from hclust
  tree
}
