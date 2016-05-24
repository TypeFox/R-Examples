SubgraphSplit <-
function(Vertices,Edges){
  # Split subgraphs using walktrap.community algorithm
  #
  # input
  #    Vertices: Vertices for network construction
  #    Edges: Edges for network construction
  #
  # output
  #    SubgroupSNPs: SNPs in the subgraphs
  #
  # Junliang Shang
  # 3.31/2014
  
  #library(igraph)
  gdf <- graph.data.frame(Edges,directed=F,Vertices)
  
  # subgraph split
  com <- walktrap.community(graph=gdf,steps=5)
  Subgroup <- split(com$names, com$membership)
  SubgroupNum <- length(Subgroup)
  SubgroupSNPs <- vector("list",SubgroupNum)
  # cat("There are ",SubgroupNum," subgroups \n\n")
  
  for (i in 1:SubgroupNum){
    # cat(i,":\n")
    AllNodes <- Subgroup[[i]]
    Whether <- sapply(X=AllNodes,CheckVerticeType)
    Whether <- sapply(1:length(AllNodes),function(x,Whether) Whether[[x]]
                      ,Whether)
    AllNodes <- AllNodes[Whether]
    if (length(AllNodes)==0){
      cat("No real Vertices, i.e., SNPs.")
    }
   # cat("\n")
    SubgroupSNPs[[i]] <- AllNodes
  }
  
  list(SubgroupSNPs=SubgroupSNPs)
}
