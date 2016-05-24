#'
#' Phylogenetic tree exploration
#'
#' Compares phylogenetic trees and maps them into a small number of dimensions for easy visualisation and identification of clusters.
#'
#' @param x an object of the class multiPhylo
#' @param method the method for summarising the tree as a vector.
#' Choose from:
#' \code{treeVec} (default) the Kendall Colijn metric vector
#' \code{RF} the Robinson Foulds metric using \code{RF.dist} from package \code{phangorn} (Note: this considers the trees as unrooted and issues a corresponding warning)
#' \code{BHV} the Billera, Holmes Vogtmann metric using \code{dist.multiPhylo} from package \code{distory}
#' others inherited from \code{distTips} in \code{adephylo}:
#' \itemize{
#' \item \code{patristic}: for each pair of tips, the sum of branch lengths on the path between them
#' \item \code{nNodes}: for each pair of tips, the number of nodes on the path between them
#' \item \code{Abouheif}: performs Abouheif's test. See Pavoine et al. (2008) and \code{adephylo}.
#' \item \code{sumDD}: sum of direct descendants of all nodes on the path, related to Abouheif's test. See \code{adephylo}.
#' }
#' @param nf the number of principal components to retain
#' @param return.tree.vectors option to also return the tree vectors. Note that this can use a lot of memory so defaults to \code{FALSE}.
#' @param ... further arguments to be passed to \code{method}.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape
#' @importFrom ade4 dudi.pco cailliez
#' @importFrom adephylo distTips
#' @importFrom distory dist.multiPhylo
#' @importFrom fields rdist
#' @importFrom phangorn RF.dist

#'
#' @examples
#'
#' ## generate list of trees
#' x <- rmtree(10, 20)
#' names(x) <- paste("tree", 1:10, sep = "")
#'
#' ## use treescape
#' res <- treescape(x, nf=3)
#' table.paint(as.matrix(res$D))
#' scatter(res$pco)
#'
#' data(woodmiceTrees)
#' woodmiceDists <- treescape(woodmiceTrees,nf=3)
#' plot(woodmiceDists$pco$li[,1],woodmiceDists$pco$li[,2])
#' woodmicedf <- woodmiceDists$pco$li
#' if(require(ggplot2)){
#' woodmiceplot <- ggplot(woodmicedf, aes(x=A1, y=A2)) # create plot
#' woodmiceplot + geom_density2d(colour="gray80") + # contour lines
#' geom_point(size=6, shape=1, colour="gray50") + # grey edges
#' geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
#' xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
#' }
#'
#' \dontrun{
#' if(require(rgl)){
#' plot3d(woodmicedf[,1], woodmicedf[,2], woodmicedf[,3], type="s", size=1.5,
#' col="navy", alpha=0.5, xlab="", ylab="", zlab="")
#' }
#' }
#'
#'
#' @export
treescape <- function(x, method="treeVec", nf=NULL, return.tree.vectors=FALSE, ...){
    ## CHECKS ##
    if(!inherits(x, "multiPhylo")) stop("x should be a multiphylo object")
    num_trees <- length(x) # number of trees
    ## fix potential bug with input of two trees
    if(num_trees<3) {
      stop("treescape expects at least three trees. The function treeDist is suitable for comparing two trees.")
    }

    # check for user supplying invalid options (these gave unhelpful error messages before)
    dots <- list(...)
    if(!is.null(dots$return.lambda.function)) stop("return.lambda.function is not compatible with treescape. Consider using multiDist instead.")
    if(!is.null(dots$save.memory)) stop("save.memory is not compatible with treescape. Consider using multiDist instead.")

    # make name labels well defined
    if(is.null(names(x))) names(x) <- 1:num_trees
    else if(length(unique(names(x)))!=num_trees){
      warning("duplicates detected in tree labels - using generic names")
      names(x) <- 1:num_trees
      }
    lab <- names(x)

    # check all trees have same tip labels
    for (i in 1:num_trees) {
      if (!setequal(x[[i]]$tip.label,x[[1]]$tip.label)) {
        stop(paste0("Tree ",lab[[i]]," has different tip labels from the first tree."))
      }
    }

    ## GET DISTANCES BETWEEN TREES, according to method ##
    ## get data.frame of all summary vectors ##
    if (method=="treeVec") {
      df <- t(data.frame(lapply(x, function(e) as.vector(treeVec(e, ...)))))
      ## get pairwise Euclidean distances ##
      D <- as.dist(rdist(df))
    }
    else if(method %in% c("patristic","nNodes","Abouheif","sumDD")){
      df <- t(data.frame(lapply(x, function(e) as.vector(adephylo::distTips(e,method=method,...)))))
      ## get pairwise Euclidean distances ##
      D <- as.dist(rdist(df))
    }
    else if(method=="RF"){
      D <- RF.dist(x)
      ## make the distance Euclidean
      D <- ade4::cailliez(D, print=FALSE)
    }
    else if(method=="BHV"){
      D <- dist.multiPhylo(x)
      ## make the distance Euclidean
      D <- ade4::cailliez(D, print=FALSE)
    }

    ## restore labels
    attr(D,"Labels") <- lab

    ## perform PCoA/MDS ##
    pco <- dudi.pco(D, scannf=is.null(nf), nf=nf)


    ## BUILD RESULT AND RETURN ##
    if (return.tree.vectors==TRUE) {
    out <- list(D=D, pco=pco, vectors=df)
    }
    else {
    out <- list(D=D, pco=pco)
    }
    return(out)
} # end treescape
