#' Function to test the significance of communities within a graph
#'
#' \code{dCommSignif} is supposed to test the significance of communities within a graph. For a community of the graph, it first calculates two types of degrees for each node: degrees based on parters only within the community itself, and the degrees based on its parters NOT in the community but in the graph. Then, it performs two-sample Wilcoxon tests on these two types of degrees to produce the signficance level (p-value)
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param comm an object of class "communities". Details on this class can be found at \url{http://igraph.org/r/doc/communities.html}
#' @return 
#' \itemize{
#'  \item{\code{significance}: a vector of p-values (significance)}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dCommSignif}}
#' @include dCommSignif.r
#' @examples
#' # 1) generate an vector consisting of random values from beta distribution
#' x <- rbeta(1000, shape1=0.5, shape2=1)
#'
#' # 2) fit a p-value distribution under beta-uniform mixture model
#' fit <- dBUMfit(x, ntry=1, hist.bum=FALSE, contour.bum=FALSE)
#'
#' # 3) calculate the scores according to the fitted BUM and fdr=0.01
#' # using "pdf" method
#' scores <- dBUMscore(fit, method="pdf", fdr=0.05, scatter.bum=FALSE)
#' names(scores) <- as.character(1:length(scores))
#'
#' # 4) generate a random graph according to the ER model
#' g <- erdos.renyi.game(1000, 1/100)
#'
#' # 5) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#'
#' # 6) find the module with the maximum score
#' module <- dNetFind(subg, scores)
#'
#' # 7) find the module and test its signficance
#' comm <- walktrap.community(module, modularity=TRUE)
#' significance <- dCommSignif(module, comm)

dCommSignif <- function(g, comm)
{

    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    if(is.null(V(ig)$name)){
        V(ig)$name <- as.character(V(ig))
    }
    
    if (class(comm) != "communities"){
        stop("The function must apply to 'communities' object.\n")
    }
    
    
    # a function to test community significance
    community.significance.test <- function(g, vids, ...) {
        subg <- igraph::induced.subgraph(g, vids)
        within.degrees <- igraph::degree(subg)
        cross.degrees <- igraph::degree(g, vids) - within.degrees
        stats::wilcox.test(within.degrees, cross.degrees, ...)
    }
    
    significance <- sapply(1:length(comm), function(x) {
        tmp <- suppressWarnings(community.significance.test(ig, vids=V(ig)$name[comm$membership==x]))
        signif(tmp$p.value, digits=3)
    })

    return(significance)
}