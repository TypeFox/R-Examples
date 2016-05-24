#'@title Caculating the strong connected components (SCC) of a network

#'@description This function utilizes Kosaraju's algorithm to caculate the
#'strong connetected components descomposition of a given network 

#'@param g, a igraph object to be caculated 
#'@references \emph{AV Aho, JE Hopcroft, JD Ullman: The design and analysis of computer algorithms, 1974}
#'@export
#'@return a list which length is equal to the number of SCCs, each element represents a Scc
#'@seealso \code{\link{getSeedSets}}
#'@examples
#'\dontrun{
#'metabolic.data <- getOrgMetabolicData("buc")
#'## metabolic network reconstruction
#'net <- reconstructGsMN(metabolic.data)
#'scc <- KosarajuSCC(net)
#'} 

KosarajuSCC <- function(g){  
  #Comparing d with S.node, get the overlapping ones and remove the NA
  getnotNA <- function(d,S.node){
    k <- length(which(d>0))   
    d.useful <- NULL
    d.useful <- subset(d[1:k],match(d[1:k],S.node)>0)
    return(d.useful)
  }
  #-----------------------------------------------------------##
  g.adjacency <- get.adjacency(g)
  n <- nrow(g.adjacency)
  List.out <- vector(mode="list")
  List.node <- seq(0,0,length.out=n)
  S.node <- seq(from=1,to=n)
##------------------declaration of variables used-------------##
  
  #g <- graph.adjacency(g.adjacency)
  #Mt <- t(M)
  gt <- t(g.adjacency) %>%
    graph.adjacency
  Source.node <- S.node[round(runif(1)*n)]
  Count <- 0
##------------------generate a stack of nodes in graph---------##
##------------------give the topology rank of each nodes-------##
  while (any(S.node != 0)){
    d <- graph.dfs(g,root=Source.node,neimode="out",unreachable=FALSE)
    d.useful <- getnotNA(d$order,S.node)
    k <- length(d.useful)
    for(i in 1:k){
      List.node[n-(Count + i-1)] <- d.useful[k-i+1]
    }
    Count <- Count+k
    S.node[d.useful] <- 0
    if(all(S.node==0)) break
    Source.node <- S.node[S.node!=0][1]
  }  
##------Reverse the directions of all arcs to obtain the transport graph---##
##-------------------------------------------------------------------------##
  Count <- 0
  while (any(List.node !=0)){
    Source.node <- List.node[which(List.node!=0)[1]]
    d <- graph.dfs(gt,root=Source.node,neimode="out",unreachable=FALSE)
    d.useful <- getnotNA(d$order,List.node)
    if(length(d)==0)
      break
    Count <- Count+1
##----Record this SCC and remove all these nodes from the graph and the stack---##
    List.out[[Count]] <- d.useful
    List.node[match(d.useful,List.node)] <- 0
  }
  List.out
}