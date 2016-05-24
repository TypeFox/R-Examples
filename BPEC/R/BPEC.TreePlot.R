BPEC.TreePlot <- function(MCMCout,colorcode=c(7,5,6,3,2,8,4,9,10))
{
  writeLines("Creating clustered tree plot...")
  count = MCMCout$countR
  root = which.max(MCMCout$RootProbsR)

  levels = MCMCout$levelsR
  datsiz = MCMCout$NoSamplesR
  clado = MCMCout$cladoR
  clusterprobs = MCMCout$ClusterProbsR
  SeqLabels =MCMCout$SeqsFileR[MCMCout$SeqLabelsR]

  if(length(SeqLabels)<count)
   {
     SeqLabels=c(SeqLabels,rep(0,count-length(SeqLabels)))
   }
  
  plotheight=sqrt(2)
 
  maxmig=dim(clusterprobs)[2]-1
  for(j in 1:count)
    {
      for(i in 1:(maxmig))
        {
          for(k in (i+1):(maxmig+1))
            {
              clusterprobs[j,i]=clusterprobs[j,i]+clusterprobs[j,k]    
            }
        }
    }
  
  for(i in 1:count)
    {
      levels[i]=levels[i]+1
    }

  NodeCoords=array(0,2*count)
  NodeOrder=array(0,count)
  newNodeOrder=array(0,count)
 
  MaxLevel=max(levels)
  plot(0,0,xlim=c(-0.01,1.01),ylim=c(-sqrt(2)*1.01,0.01),type="n",axes=FALSE,xlab="",ylab="")
  NodeCoords[2*(root-1)+1]=1/2
  NodeCoords[2*(root-1)+2]=0
  NodeOrder[1]=root
  newNodeOrder[1]=root
  descendantcounter=2
  descendants=0

  for(i in 1:MaxLevel)
    {
      prevdescendants=descendants
      
      descendants=1
      NodeOrder=newNodeOrder
      for(j in 1:count)
        {
          if(levels[j]==i)
            {
              descendants=descendants+1
            }
        }
      prevord=descendantcounter-1
      descendantcounter=1
      previousone=-1
      
      if(prevord>0)
        {
          for(j in 1:prevord)
            {           
              for(l in 1:count)
                {
                  if((clado[(NodeOrder[j]-1)*count+l]==1||clado[(l-1)*count+NodeOrder[j]]==1)&&levels[l]==i)
                    {
                      NodeCoords[(l-1)*2+1]=descendantcounter/descendants
                      NodeCoords[(l-1)*2+2]=-levels[l]/MaxLevel*plotheight
                      if((NodeCoords[(NodeOrder[j]-1)*2+1]<(descendantcounter+1)/descendants&&(NodeCoords[(NodeOrder[j]-1)*2+1]>NodeCoords[(previousone-1)*2+1]||previousone==-1))||descendants==1)
                        {
                          NodeCoords[(l-1)*2+1]=NodeCoords[(NodeOrder[j]-1)*2+1]
                        }

                      thickness = MCMCout$EdgeTotalProbR[NodeOrder[j],l] * 1.5
                      lines(c(NodeCoords[(l-1)*2+1],NodeCoords[(NodeOrder[j]-1)*2+1]),c(NodeCoords[(l-1)*2+2],NodeCoords[(NodeOrder[j]-1)*2+2]),lwd=thickness)
                    
                      previousone=l
                      newNodeOrder[descendantcounter]=l
                      descendantcounter=descendantcounter+1
                      if(i==MaxLevel)
                        {
                          if(datsiz[l]>0)
                            {
                              for(w in 1:(maxmig+1))
                                {
                                  BPEC.circles(x=NodeCoords[(l-1)*2+1],y=NodeCoords[(l-1)*2+2],radius=0.02+0.0005*datsiz[l],col=colorcode[w],angle=360*(clusterprobs[l,w]))
                                }
                              text(NodeCoords[(l-1)*2+1],NodeCoords[(l-1)*2+2],SeqLabels[l],cex=1)
                            }
                          if(datsiz[l]==0)
                            {
                              BPEC.circles(x=NodeCoords[(l-1)*2+1],y=NodeCoords[(l-1)*2+2],radius=0.003,col=1)
                            }
                        }
                      
                      next
                    }
                }
#####
              if(datsiz[NodeOrder[j]]>0)
                {                 
                  for(w in 1:(maxmig+1))
                    {                   
                      BPEC.circles(x=NodeCoords[(NodeOrder[j]-1)*2+1],y=NodeCoords[(NodeOrder[j]-1)*2+2],radius=0.02+0.0005*datsiz[NodeOrder[j]],col=colorcode[w],angle=360*(clusterprobs[NodeOrder[j],w]))                   
                    }
               
                  text(NodeCoords[(NodeOrder[j]-1)*2+1],NodeCoords[(NodeOrder[j]-1)*2+2],SeqLabels[NodeOrder[j]],cex=1)
                }
              if(datsiz[NodeOrder[j]]==0)
                {
                  BPEC.circles(x=NodeCoords[(NodeOrder[j]-1)*2+1],y=NodeCoords[(NodeOrder[j]-1)*2+2],radius=0.003,col=1)
                }
#####
            }
        }
    }

  Output=list()
#  writeLines("Creating GoogleEarth Tree plot...")
  TreeEdges = MCMCout$TreeEdges
  clustprob = MCMCout$ClusterProbsR
  count = MCMCout$countR
  dims = dim(MCMCout$SampleMeansR)[1]
####################################################################
                                        # required libraries igraph, R2G2, ape
####################################################################
                                        # include network.to.newick.r - function
####################################################################
                                        #source("network.to.newick.mod.r")
                                        #source("~/Desktop/BPEC-Rsources/network.to.newick_igraph.r")
####################################################################
                                        # load needed library "igraph"
                                        # library("igraph")
                                        # make graph from edgelist
                                        #TreeEdgesOut = data.matrix(TreeEdges[,1:2])
  TreeEdgesOut = data.matrix(MCMCout$TreeEdges[,1:2])
  
  dimnames(TreeEdgesOut) = NULL
  GraphEdges = graph.edgelist(TreeEdgesOut, directed=TRUE)
                                        # name vertex sequence
                                        #V(GraphEdges)$name  = paste("h", sep="", V(GraphEdges))
  V(GraphEdges)$name  = paste(V(GraphEdges))
                                        # remove un-connected vertices
  GraphEdgesSub = subgraph.edges(graph=GraphEdges, eids=1:length(E(GraphEdges)), delete.vertices = TRUE)
                                        #GraphEdgesSub
  roundInt = 1000 * round(clustprob, 3)
    rowMat = split(roundInt, row(roundInt))
    attributes(rowMat) = NULL
    VertShape = ifelse(!is.na(roundInt[c(1:count)])==T, "pie", "none")
  
                                        #call the pdf writer
                                        # pdf(file="HaplotypeSubgraphNetwork.pdf", paper="a4r", width = 0, height = 0)
    par(mai=c(0,0,1,0))
    par(mfrow=c(1,2)) 
    set.seed(12345)
    
    plot.igraph(GraphEdgesSub, layout=layout.kamada.kawai, 
                vertex.shape=VertShape, vertex.pie=rowMat,
                vertex.pie.color=list(colorcode), vertex.pie.lty=0,
                vertex.size=6, vertex.label=V(GraphEdgesSub)$name,
                vertex.label.cex=0.6, vertex.label.font=2,vertex.label.dist=0,
                edge.color="black", edge.arrow.size = 0.05)
    ##run the plot
                                        #  dev.off() #close the device   
    Output$GraphEdgesSub = GraphEdgesSub

                                  #TreeEdgesOut = data.matrix(TreeEdges[,1:2])
    TreeEdgesOut = data.matrix(MCMCout$TreeEdges[,1:2])
    
    dimnames(TreeEdgesOut) = NULL
    GraphEdges = graph.edgelist(TreeEdgesOut, directed=TRUE)
                                        # name vertex sequence
                                        #V(GraphEdges)$name  = paste("h", sep="", V(GraphEdges))
    V(GraphEdges)$name  = paste(V(GraphEdges))
                                        # remove un-connected vertices
    GraphEdgesSub = subgraph.edges(graph=GraphEdges, eids=1:length(E(GraphEdges)), delete.vertices = TRUE)
                                        #GraphEdgesSub
    
                                        # preparation for haplotype-graph plotting
    clustprob[clustprob %in% NaN] = NA
    
                                        # make proportions (rounded integers that sum up to 1000)
                                        #roundint = round(clusterprobs * datsiz[1:count])
    roundInt = 1000 * round(clustprob, 3)
    rowMat = split(roundInt, row(roundInt))
    attributes(rowMat) = NULL
    
                                        # create newick string without lengths
                                        #GraphEdges.nwk = network.to.newick(GraphEdgesSub)
                                        # or
    GraphEdges.nwk = explode(GraphEdgesSub)
                                        # string manipulation
    GraphEdges.nwk = paste("(",strsplit(GraphEdges.nwk,"\\;"),");", sep="")
    
    ## input newick string to create a tree
    
    GraphEdgesTree = read.newick(text=GraphEdges.nwk)
                                        # remove singletons
    GraphEdgesTree = collapse.singles(GraphEdgesTree)

                                        #   pdf(file="HaplotreeNoSingles.pdf", paper="a4r", width = 0, height = 0)
    plot(GraphEdgesTree, type = "c", 
         direction = "l", adj = 1,
         show.node.label = T, label.offset = 0.5, srt = 20,
         font = 3, cex = 0.4, 
         root.edge = T)
                                        #   dev.off()
    Output$GraphEdgesTree = GraphEdgesTree
  return(Output)
}

