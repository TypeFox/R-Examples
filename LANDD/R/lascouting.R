#' Find the liquid association scouting genes.
#' 
#' \code{lascouting()} is used to find the liquid association scouting gene 
#' @useDynLib LANDD
#' @param network.graph An igraph object representing the gene network.
#' @param express.matrix A matrix representing the expression matrix for the genes in gene 
#' network. Row names are the gene ids in gene network.
#' @param k Integer giving the order of the ego-network.
#' @param n.cores Number of cores used for parallel computing.
#' @return A logical matrix representing the LA-scouting genes for each gene. Rows represent 
#' the ego gene id and columns represents the LA-scouting genes.
#' @examples \dontrun{laresult <- lascouting(g,m,k=2,n.cores=4)}
#' @export
#' @import igraph
#' @importFrom Matrix Matrix
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom fdrtool fdrtool
#' @importFrom Matrix Matrix
#' @importFrom Rcpp evalCpp
#' 
#' @importFrom modeest mlv
#' 
#' 
lascouting <- function(network.graph, express.matrix, k = 2, n.cores = 4) {
  network.node <- V(network.graph)$name
  matrix.node <- row.names(express.matrix)
  if (is.null(network.node) || is.null(matrix.node)) 
    stop("node name can't be null")
  if (!identical(intersect(network.node, matrix.node), union(network.node, matrix.node))) {
    common.node <- getCommonNode(network.graph, express.matrix)
    network.graph <- cleanGraph(network.graph, common.node)
    express.matrix <- cleanMatrix(express.matrix, common.node)
  } else {
    common.node <- network.node
  }
  size <- length(common.node)
  
  express.matrix <- normalizeInputMatrix(express.matrix)
  
  if (k != 1) {
    graph.connected <- connect.neighborhood(network.graph, k)
    connected.list <- as.matrix(get.edgelist(graph.connected))
  } else {
    connected.list <- as.matrix(get.edgelist(network.graph))
  }
  row.size <- nrow(connected.list)
  express.matrix.t <- t(express.matrix)/ncol(express.matrix)
  
  cl <- makeCluster(n.cores, outfile = "")
  registerDoParallel(cl)
  
  result <- foreach(i = 1:row.size) %dopar% {
    print(i)
    if (connected.list[i, 1] != connected.list[i, 2]){
      xy <- express.matrix[connected.list[i, 1], ] * express.matrix[connected.list[i, 2], ]
      la.vector <- c(xy %*% express.matrix.t)
      la.vector <- la.vector - mlv(la.vector,method="shorth")$M
      lfdr <- fdrtool(la.vector, verbose = FALSE, plot = FALSE)$lfdr
      return(rownames(express.matrix)[which(lfdr < 0.2)])
    }
    else {
      return(NULL)
    }
  }
  stopCluster(cl)
  matrix.rowname = rownames(express.matrix)
  node.z <- Matrix(0, nrow = size, ncol = size, dimnames = list(matrix.rowname, matrix.rowname))
  for (i in 1:row.size) {
    if (length(result[[i]]) != 0) {
      x <- connected.list[i, 1]
      y <- connected.list[i, 2]
      node.z[x, c(result[[i]])] <- 1
      node.z[y, c(result[[i]])] <- 1
    }
  }
  diag(node.z) <- 0
  return(node.z)
  
}