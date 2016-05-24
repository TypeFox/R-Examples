

#' Hierarchical Clustering with Prototypes: Minimax Linkage.
#' 
#' Functions to perform minimax linkage hierarchical clustering and to cut such
#' trees to return clusterings with prototypes.
#' 
#' \tabular{ll}{ Package: \tab protoclust\cr Type: \tab Package\cr Version:
#' \tab 1.0\cr Date: \tab 2011-06-21\cr License: \tab GPL-2\cr LazyLoad: \tab
#' yes\cr }
#' 
#' @name protoclust-package
#' @docType package
#' @useDynLib protoclust
#' @author Jacob Bien and Rob Tibshirani
#' 
#' Maintainer: Jacob Bien <jbien@@cornell.edu>
#' @seealso \code{\link{protoclust}}, \code{\link{protocut}},
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
NULL



