#' Find weights based on kernel density on the graph.
#'  
#' There are three common ways to invoke \code{graph.kd}:
#' \itemize{
#'   \item \code{graph.kd(relate_matrix, graph, smoothing.normalize=c('one'))}
#'   \item \code{graph.kd(relate_matrix, graph, smoothing.normalize=c('squareM'))}
#'   \item \code{graph.kd(relate_matrix, graph, smoothing.normalize=c('none'))}
#'   }
#' The first method is used when the total weight of all genes z is set to 'one'.
#' In this way, those genes surrounded by more genes z will not take advantages over those surrounded by fewer genes.
#' In contrast, the second method takes the number of genes around into consideration, the result of the first method will
#' multiply the square of the number of genes around.
#' The third method does not normalize the data. Thus genes with more neighbors are more likely to receive higher weights. 
#' 
#' @param relate.matrix The matrix returned by lascouting.
#' @param network.graph The igraph object representing the gene network.
#' @param smoothing.normalize Different ways to normalize the result.
#' @return A matrix representing the weights calculated using kernel density for each gene. Each row is an ego gene, columns
#' are the weights of potential scouting genes for the gene. 
#' @examples \dontrun{
#' relate.matrix <- lascouting(g,m,k=2,n.cores=4) 
#' graph.kd(relate.matrix,g,smoothing.normalize = "one")}
#' @export
#' @import igraph
#' @importFrom stats dnorm

graph.kd <- function(relate.matrix, network.graph, smoothing.normalize = c("one", "squareM", "none")) {
  smoothing.normalize <- match.arg(smoothing.normalize)
  
  network.node <- V(network.graph)$name
  matrix.node <- row.names(relate.matrix)
  if (!identical(intersect(network.node, matrix.node), union(network.node, matrix.node))) {
    common.node <- getCommonNode(network.graph, relate.matrix)
    network.graph <- cleanGraph(network.graph, common.node)
  }
  #weight = c(dnorm(0),dnorm(1),dnorm(2))
  weight0 <- dnorm(0)
  weight1 <- dnorm(1)
  weight2 <- dnorm(2)
  
  size <- nrow(relate.matrix)
  
  
  relate.matrix <- relate.matrix[order(rownames(relate.matrix)), ]
  relate.matrix <- relate.matrix[, order(colnames(relate.matrix))]
  
  adjacency1 <- get.adjacency(network.graph, type = "both")
  adjacency2 <- get.adjacency(connect.neighborhood(network.graph, 2), type = "both") - adjacency1
  
  adjacency1 <- adjacency1[order(rownames(adjacency1)), ]
  adjacency1 <- adjacency1[, order(colnames(adjacency1))]
  adjacency2 <- adjacency2[order(rownames(adjacency2)), ]
  adjacency2 <- adjacency2[, order(colnames(adjacency2))]
  
  temp <- diag(size) * weight0
  weight.matrix <- temp + adjacency1 * weight1 + adjacency2 * weight2
  
  if (smoothing.normalize == "one") {
    rsum <- rowSums(as.matrix(weight.matrix))
    nmatrix <- diag(1/rsum)
    colnames(nmatrix) <- rownames(weight.matrix)
    weight.matrix <- nmatrix %*% weight.matrix
  } else if (smoothing.normalize == "squareM") {
    temp.weight.matrix <- as.matrix(weight.matrix)
    rsum <- rowSums(temp.weight.matrix)
    m <- rowSums(temp.weight.matrix != 0)
    nmatrix <- diag(sqrt(m)/rsum)
    colnames(nmatrix) <- rownames(weight.matrix)
    weight.matrix <- nmatrix %*% weight.matrix
    
  }
  
  result <- relate.matrix %*% weight.matrix
  return(result)
  
}