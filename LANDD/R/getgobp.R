#' @title Record genes W
#' 
#' @description \code{get.W()} creates a table to record Gene Ontology Biological Process mapping results.  Every gene x takes a row.
#' 
#' @details \code{get.W()} generates a result file of ego gene X, genes within k steps of X, 
#' the liquid association scouting genes of x and genes W.Every gene x takes a row in the table.
#' @param graph The graph of gene network.
#' @param laresult The result of lascouting which finds the liquid association scouting genes.
#' @param z.matrix A matrix representing gene Z (selected scouting genes). Row names are the gene id in gene network.
#' @param k An Integer giving the order of the network.
#' @param cutoff The threshold to find LA scouting genes. 
#' @return A table records the intermediate result of Gene Ontology Biological Process which contains ego gene X, genes within k steps of X, 
#' the liquid association scouting genes of x and genes W.Each x occupies a row. 
#' @examples \dontrun{ 
#' laresult <- lascouting(g,m,k=2,n.cores=4)  
#' get.W(g,laresult,z,cutoff=0.8,k=2)}
#' @export
#' @importFrom stats median
#' 
get.W <- function(graph, laresult, z.matrix, cutoff, k=2) {
  xlist = row.names(z.matrix)
  LANDDdata <- setRefClass("LANDDdata",
                         fields = list(x = "character", y = "character"))
  LANDDList<-c()
  for(x in xlist){
    if(sum(laresult[x,])!=0 && length(z.matrix[x,][z.matrix[x,]>cutoff])!=0){
      print(x)
      y <- V(graph)$name[unlist(igraph::neighborhood(graph, k, nodes = x))]
      y <- paste(y[y!=x],collapse = " ")
      
      z = paste(names(laresult[x,][laresult[x,]==1]),collapse = " ")
      w = z.matrix[x,]
      w = w[w>cutoff]
      w = paste(paste(names(w),w,sep=":"),collapse = " ")
      LANDDList <- rbind(LANDDList,cbind(x,y,z,w))
    }
  }
  colnames(LANDDList)<-c("x","y","z","w:w_value")
  return(LANDDList)

  
  
} 
  
  








#' Create a table to record Gene Ontology Biological Process mapping results.  Every gene W's community takes a row.
#' 
#' \code{getgobp.community()} generates a result file of ego gene X,  significant GO terms of X, significant GO terms 
#' of genes within k steps of X, gene W, significant GO terms  of W,
# ' the similarity of gene W and genes within k steps of gene X, the average distance between gene X and gene
# W. ' A gene X may correspond with several W communities. Thus one community takes a row in the table.
#' @param graph The graph of gene network.
#' @param z.matrix A matrix representing gene Z (selected scouting genes). Row names are the gene id in gene network.
#' @param k An Integer giving the order of the network.
#' @param n.cores The number of cores used for parallel computing.
#' @param cutoff The threshold to find LA scouting genes.
#' @param community Boolean. Whether compute the community of genes W or not.
#' @param community.min Integer. The minimum number of genes numbers in a community. 
#' @param term.limit The maximum number of GO terms to list in a row of the table.
#' @return A table containing the IDs of scouting center genes W, over-represented GO terms by
#'  W, semantic similarity on the Gene Ontology system between the X ego network and all 
#'  scouting center genes, average graph distance between gene X and W. W are grouped by 
#'  network community. Each W community occupies a row. 
#' @examples \dontrun{
#' g <- graph.data.frame(as.matrix(read.table("HumanBinaryHQ_HINT.txt"))) 
#' getgobp(g,z,k=2,n.cores=4,cutoff=1,community=TRUE,community.min=5,term.limit = NA)}
#' @export
#' @import igraph
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @import GOstats
#' @import GOSemSim
getgobp <- function(graph, z.matrix, k = 2, n.cores = 4, cutoff = 1,community = TRUE, community.min = 5, term.limit = NA) {
  resulttable = NULL
  cl <- makeCluster(n.cores, outfile = "")
  all.entrez <- colnames(z.matrix)
  registerDoParallel(cl)
  if(community) {
    community <- apply(z.matrix, 1, getCommunity, graph, cutoff, community.min)
    community <- community[!sapply(community, is.null)]
    
    cat("loop begin\n")
    
    resulttable <- foreach(i = 1:length(names(community)), .combine = "rbind") %dopar% {
      x <- names(community)[i]
      wc <- community[[x]]
      member <- membership(wc)
      
      community_index <- names(sizes(wc)[sizes(wc) > community.min])
      sel.entrez <- x
      xgo <- getGO(sel.entrez, all.entrez)
      
      if (is.null(xgo) || is.na(xgo$Pvalue) || length(xgo$Term) == 0) {
        return(NULL)
      } else {
        if (!is.na(term.limit)) {
          xgo <- xgo[1:term.limit, ]
        }
        xgo <- paste(xgo$Term, signif(xgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
      }
      xk <- V(graph)[unlist(igraph::neighborhood(graph, k, nodes = x))]$name
      
      sel.entrez <- xk
      xkgo <- getGO(sel.entrez, all.entrez)
      
      
      if (is.null(xkgo) || is.na(xkgo$Pvalue) || length(xkgo$Term) == 0) {
        return(NULL)
      } else {
        if (!is.na(term.limit)) {
          xkgo <- xkgo[1:term.limit, ]
        }
        xkgo <- paste(xkgo$Term, signif(xkgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
      }
      
      w.result <- do.call("rbind", lapply(community_index, get.W.GO, member, xk, x, graph, all.entrez, term.limit))
      if (is.null(w.result)) {
        return(NULL)
      } else {
        print(x)
        return(rbind(resulttable,cbind(x, xgo, xkgo, w.result)))
      }
    }
    stopCluster(cl)
    return(resulttable)
  }
  else {

    resulttable <- foreach(i = 1:nrow(z.matrix), .combine = "rbind") %dopar% {
      
      x = rownames(z.matrix)[i]
      w <- names(z.matrix[i,z.matrix[i,]>cutoff])
      if(length(w)==0) {return(NULL)}
      xgo <- getGO(x, all.entrez)
      if (is.null(xgo) || is.na(xgo$Pvalue) || length(xgo$Term) == 0) {
        return(NULL)
      } else {
        if (!is.na(term.limit)) {
          xgo <- xgo[1:term.limit, ]
        }
        xgo <- paste(xgo$Term, signif(xgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
      }
      cat("2:",x,"\n")
      
      xk <- V(graph)[unlist(igraph::neighborhood(graph, k, nodes = x))]$name
      xkgo <- getGO(xk, all.entrez)
      
      if (is.null(xkgo) || is.na(xkgo$Pvalue) || length(xkgo$Term) == 0) {
        return(NULL)
      } else {
        if (!is.na(term.limit)) {
          xkgo <- xkgo[1:term.limit, ]
        }
        xkgo <- paste(xkgo$Term, signif(xkgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
      }
      

      wgo <- getGO(w,all.entrez)
      if (is.null(wgo) || is.na(wgo$Pvalue) || length(wgo$Term) == 0) {
        return(NULL)
      } else {
        if (!is.na(term.limit)) {
          wgo <- wgo[1:term.limit, ]
        }
        wgo <- paste(wgo$Term, signif(wgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
      }
      
      

      xk.w.semantic.similarity <- clusterSim(c(w), c(xk), combine = "avg")
      x.w.avg.distance <- mean(igraph::shortest.paths(graph, v = w, to = x))
      if(is.infinite(x.w.avg.distance))
      {
        x.w.avg.distance <- median(igraph::shortest.paths(graph, v = w, to = x))
      }
      
      w <- paste(w, collapse = " ")
      print(x)
      return(cbind(x, xgo, xkgo, w, wgo, xk.w.semantic.similarity, x.w.avg.distance))
      
      
    }
    stopCluster(cl)
    return(resulttable)
  }
  
}
globalVariables('i')




get.W.GO <- function(ci, member, xk, x, graph, all.entrez, term.limit) {
  
  
  w <- names(member[member == ci])
  
  sel.entrez <- w
  wgo <- getGO(sel.entrez, all.entrez)
  if (is.null(wgo) || is.na(wgo$Pvalue) || length(wgo$Term) == 0) {
    return(NULL)
  } else {
    if (!is.na(term.limit)) {
      wgo <- wgo[1:term.limit, ]
    }
    wgo <- paste(wgo$Term, signif(wgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
    
  }
  
  
  xk.w.semantic.similarity <- clusterSim(c(w), c(xk), combine = "avg")
  x.w.avg.distance <- mean(igraph::shortest.paths(graph, v = w, to = x))
  if(is.infinite(x.w.avg.distance))
  {
    x.w.avg.distance <- median(igraph::shortest.paths(graph, v = w, to = x))
  }
  w <- paste(w, collapse = " ")
  
  return(data.frame(w, wgo, xk.w.semantic.similarity, x.w.avg.distance))
}




