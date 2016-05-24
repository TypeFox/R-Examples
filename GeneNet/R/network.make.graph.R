### network.make.graph  (2015-08-02)
###
###   Construct "graph" object from given network
###
### Copyright 2003-15 Juliane Schaefer and Korbinian Strimmer
###
###
### This file is part of the `GeneNet' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# requires the installation of the "graph" library

# generate a graph object from and edge list
# (such as obtained from network.test.edges) 
network.make.graph = function(edge.list, node.labels, drop.singles=FALSE)
{
    V = unique(node.labels)  
    if ( length(V) != length(node.labels) )
    {
       stop("Duplicate node labels encountered. Node labels must be unique!")
    }   
    V = as.character(V)
     
    
    # create empty graph with no edges
    gX = graph::graphNEL(nodes=V, edgemode="directed") 

    # check whether we might have directed edges
    if (ncol(edge.list) == 11)
      undirected=FALSE
    else
      undirected=TRUE

    # add edges and edge weights (correlations)
    sourceNode=NULL
    destNode=NULL
    edgeCorr=NULL
    for(i in 1:nrow(edge.list))   
    {
        rc = round(edge.list[i,1], digits=5)
        if (undirected)
        {
          # add edge twice in both directions
          sourceNode=c(V[edge.list[i,2]], V[edge.list[i,3]],  sourceNode)
          destNode  =c(V[edge.list[i,3]], V[edge.list[i,2]],  destNode)
          edgeCorr=c( rc, rc, edgeCorr )
        }
        else
        {
            if(edge.list[i,11]=="1to2" | edge.list[i,11]=="undirected")
            {
                # add single edge in direction 1->2
                sourceNode=c(V[edge.list[i,2]],  sourceNode)
                destNode  =c(V[edge.list[i,3]],  destNode)
                edgeCorr=c( rc, edgeCorr )   
            }

            if (edge.list[i,11]=="2to1" | edge.list[i,11]=="undirected")
            {
                 # add single edge in direction 2->1
                sourceNode=c(V[edge.list[i,3]],  sourceNode)
                destNode  =c(V[edge.list[i,2]],  destNode)
                edgeCorr=c( rc, edgeCorr )   
            }           
        }
    }
    gX = graph::addEdge(sourceNode, destNode, gX, edgeCorr )


    if(drop.singles) # remove unconnected nodes
    {
      # nodes with degree > 0
      nd = graph::nodes(gX)[ node.degree(gX) > 0 ]
    
      gX = graph::subGraph(nd, gX)
    }
  
    return(gX)
}

# number of edges connected to a node (bi-directional edges are counted only once)
node.degree = function(gr)
{
  return( graph::degree(graph::ugraph(gr)))
}

# show number of nodes
num.nodes = function(gr)
{
  return( length(graph::nodes(gr)) )
}


# print edge weights and edge directions
edge.info = function(gr)
{
    em = graph::edgeMatrix(gr, duplicates=FALSE) 
    ned = dim(em)[2]
    ds = rep("forward", ned)

    # for (partially) directed graphs, remove one of the edges for undirected edges
    duplicatedEdge = rep(FALSE, ned)
    
    for (i in 1:ned)
    {
      # check for duplication only if not already recognized as duplicated edge
      if(!duplicatedEdge[i]) 
      {
        reverseEdge = c(em[2,i],em[1,i])
        dupidx = which( apply(reverseEdge == em , 2, function(x) x[1] && x[2] ))
        if (length(dupidx) > 0)
        {
          duplicatedEdge[dupidx] = TRUE
          ds[i] = "none"
        }
      }
    }

    em2 = em[, !duplicatedEdge]
    weight = graph::eWV(gr, em2, useNNames = TRUE, sep="~") # edge weights

    dir = ds[ !duplicatedEdge]   # edge direction
    names(dir) = names(weight)
      
    return( list( weight=weight, dir=dir) )
}

