getCommonNode <- function(network.graph, matrix) {
  network.node <- V(network.graph)$name
  matrix.node <- row.names(matrix)
  common.node <- intersect(network.node, matrix.node)
  return(common.node)
}

cleanGraph <- function(network.graph, remain) {
  network.node <- V(network.graph)$name
  delete <- setdiff(network.node, remain)
  network.graph <- igraph::delete.vertices(network.graph, delete)
  return(network.graph)
}

cleanMatrix <- function(express.matrix, remain) {
  express.matrix <- express.matrix[rownames(express.matrix) %in% remain, ]
  return(express.matrix)
}

lascore <- function(x, y, z) {
  sum(x * y * z)/length(x)
}

getSubNet <- function(graph, x, order) {
  return(neighborhood(graph, order, x))
}


getCommunity <- function(z, g, cutoff, community.min) {
  zcutoff <- names(z[z > cutoff])
  if (length(zcutoff) == 0) {
    return(NULL)
  }
  subg <- induced.subgraph(graph = g, vids = zcutoff)
  wc <- igraph::walktrap.community(subg)
  
  member <- membership(wc)
  w <- names(member[member == which.max(sizes(wc))])
  if (length(w) >= community.min) {
    return(wc)
  } else {
    return(NULL)
  }
}


getGO <- function(sel.entrez, all.entrez) {
  params <- new("GOHyperGParams", geneIds = sel.entrez, universeGeneIds = all.entrez, ontology = "BP", pvalueCutoff = 0.01, 
                conditional = F, testDirection = "over", annotation = "hgu133plus2.db")
  Over.pres <- tryCatch({
    Over.pres <- hyperGTest(params)
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(Over.pres)) {
    return(NULL)
  }
  
  summary <- getGeneric("summary")
  ov <- summary(Over.pres)
  return(ov[ov$Size < 1000 & ov$Size > 5, ])
}

cutoffz <- function(z, cutoff) {
  zcutoff <- names(z[z > cutoff])
  if (length(zcutoff) == 0) {
    return(NULL)
  }
  return(zcutoff)
}

#' Create a table to record the distance between gene x and gene w.
#' 
#' \code{xw.distance()} generates a table contains distance between all genes x and their correspongding genes w. 
#' @param graph The graph of the gene network.
#' @param z.matrix A matrix representing gene Z. Row names are the gene id in gene network. 
#' @param cutoff A number used to find LA scouting gene z. 
#' @param n.cores Core number used for parallel computing.
#' @return a table contains distance between all genes x and their correspongding genes w. 
#' @examples \dontrun{ xw.distance(g,m,cutoff=0.8,n.cores=4)}
#' @export
#' 
#' 
xw.distance <- function(graph, z.matrix, cutoff = 0.8, n.cores = 4) {
  wlist <- apply(z.matrix, 1, cutoffz, cutoff)
  wlist <- wlist[!sapply(wlist, is.null)]
  
  cl <- makeCluster(n.cores, outfile = "")
  registerDoParallel(cl)
  resulttable <- foreach(j = 1:length(names(wlist)), .combine = "rbind") %dopar% {
    x <- names(wlist)[j]
    wl <- wlist[[x]]
    distance.table <- NULL
    for (w in wl) {
      distance <- shortest.paths(graph, v = x, to = w)
      
      distance.table <- rbind(distance.table, c(x, w, distance))
    }
    return(distance.table)
  }
  
}

globalVariables('j') 

w.distance <- function(ci, member, xk) {
  w <- names(member[member == ci])
  sim <- clusterSim(c(w), c(xk), combine = "avg")
  w <- paste(as.character(w), collapse = " ")
  xk <- paste(xk, collapse = " ")
  return(data.frame(xk, w, sim))
}