traceNetworks <-
function(ARTIVAnet, edgesThreshold, parentColor = "blue", targetColor = "grey", parentgeneNames = TRUE, targetgeneNames = TRUE, layout = "fruchterman.reingold", onepage = TRUE)
{
  # Input parameters:
  # ARTIVAnet: is a data frame with the results obtained after an ARTIVA
  # analysis.
  # Threshold: is a value to select an interaction between target and
  # parent genes (based posterior probabilities).

  if(edgesThreshold<min(ARTIVAnet$edgesThreshold)){
    print(paste("WARNING : The coefficients for edges with posterior probability below", min(ARTIVAnet$edgesThreshold), "were not estimated (grey edges) in the network to be plotted."))
  }
  
  # Lists of different target genes and parent genes to be represented
  targetGeneList = as.character(unique(ARTIVAnet$Target))
  parentGeneList = as.character(unique(ARTIVAnet$Parent))
  
  # Number of genes in the graphical
  GeneNumber = length(targetGeneList) + length(parentGeneList)
  # Create a graph with all parent and target genes as nodes
  
  # To escape missing values
  ARTIVAnet$PostProb[which(is.na(ARTIVAnet$PostProb))] = 0
    
  # Create the graph with all potential interactions between genes
  PotentialEdgesList = cbind(as.character(ARTIVAnet$Parent), 
			     as.character(ARTIVAnet$Target))
  GlobalNetwork = graph.edgelist(PotentialEdgesList)

  # CoeffMean values are added to ponderate edges (positive or negative values, etc.)
  GlobalNetwork = set.edge.attribute(GlobalNetwork, "weight", 
				      value = ARTIVAnet$CoeffMean) 
  # Color for the arrows to be drawn
  EdgeColorVec = rep("black", length(ARTIVAnet$CoeffMean))
  # Induction = red and repression = green
#  colInduction="#B9121B"
#  colRepression="#96CA2D"
#  EdgeColorVec[ARTIVAnet$CoeffMean > 0] = colInduction
#  EdgeColorVec[ARTIVAnet$CoeffMean < 0] = colRepression
  E(GlobalNetwork)$color = EdgeColorVec

  # Color for the arrows to be drawn
  EdgeTypeVec = rep(1, length(ARTIVAnet$CoeffMean))
  # Induction = red and repression = green
  TypeInduction=1
  TypeRepression=2
  EdgeTypeVec[ARTIVAnet$CoeffMean > 0] = TypeInduction
  EdgeTypeVec[ARTIVAnet$CoeffMean < 0] = TypeRepression

  E(GlobalNetwork)$lty = EdgeTypeVec

#  GlobalNetwork = set.edge.attribute(GlobalNetwork, "color", 
#				      value = EdgeColorVec)

  # Keep only edges in the graph for which posterior probability is 
  # upper the threshold
  EdgesToDelete = which(ARTIVAnet$PostProb < edgesThreshold)
  SubNetwork = delete.edges(GlobalNetwork, E(GlobalNetwork)[EdgesToDelete - 1])
  # Note that indice in the vector E(g) start at 0


########################################################################
#Modified: 02/17/2012
#Last modification: 03/16/2012
#By SERVILLO David
########################################################################
#***************************************************************************************************
  # Names of the nodes represented in the graph.
	NodeLabel  = get.vertex.attribute(GlobalNetwork, name = "name")

	#Declare and initialize intervals of parent genes and target genes
	intervalParent = 0
	intervalTarget = 0

	SizeParent = 20
	SizeTarget = 10

	NodeCoord = NULL
  # Nodes coordinates are calculated according to the global structure of the graph
	if(layout == "geneLines")   #The user did not specify the layout
	{

		#The space beetwen vertices on the same line are calculated
		intervalParent = 2/(length(parentGeneList) - 1)
		intervalTarget = 2/(length(targetGeneList) - 1)
		xParent = -1
		xTarget = -1
		for(i in 1:length(NodeLabel))
		{
			if(sum(parentGeneList == NodeLabel[i]) > 0)   #Parent genes are drawn in a line at the top of the graphical interface
			{
				if(length(parentGeneList) == 1)
					NodeCoord=rbind(NodeCoord, c(0, 0.8))
				else
				{
					NodeCoord=rbind(NodeCoord, c(xParent, 0.8))
					xParent = xParent + intervalParent
				}
			}
			else if(sum(targetGeneList == NodeLabel[i]) > 0)   #Target genes are drawn in a line at the bottom of the graphical interface
			{
				if(length(targetGeneList) == 1)
				{
					NodeCoord=rbind(NodeCoord, c(0, 0.2))
				}
				else
				{
					NodeCoord=rbind(NodeCoord, c(xTarget, 0.2))
					xTarget = xTarget + intervalTarget
				}
			}
		}
	#The size of the vertices are reduced if the interval is too low
	if(intervalParent < 0.2)
	{SizeParent = 210/length(parentGeneList)}
	if(intervalTarget < 0.1)
	{SizeTarget = 200/length(targetGeneList)}

    # end of plot with geneLines()
	}


  # Nodes coordinates are calculated according to the global structure of the graph  
  if(layout == "fruchterman.reingold")
  {
	  NodeCoord  = layout.fruchterman.reingold(SubNetwork) 
  }
  if(layout == "random")
    {
      NodeCoord  = layout.random(SubNetwork) 
    }
  if(layout == "circle")
    {
      NodeCoord  = layout.circle(SubNetwork) 
    }  		
  if(layout == "kamada.kawai")
    {
      NodeCoord  = layout.kamada.kawai(SubNetwork)
    }
  if(layout == "spring")
    {
      NodeCoord  = layout.spring(SubNetwork)
    }
  if(layout == "reingold.tilford")
    {
      NodeCoord  = layout.reingold.tilford(SubNetwork)
    }
  if(layout == "lgl")
    {
      NodeCoord  = layout.lgl(SubNetwork)
    }  	
  if(layout == "graphopt")
    {
      NodeCoord  = layout.graphopt(SubNetwork)
    }  	
  if(layout == "mds")
    {
      NodeCoord  = layout.mds(SubNetwork)
    }		
  if(layout == "svd")
    {
      NodeCoord  = layout.svd(SubNetwork)
    }  	

  # Names of the nodes represented in the graph.
  NodeLabel  = get.vertex.attribute(GlobalNetwork, name = "name")

  # Vectors to store different information for the graphical representation
  NodeSize = NULL
  NodeColor = NULL
  #LabelPos = NULL
  
  # Labels of the nodes to be written
  WritenNodeLabel = NodeLabel

  #To adapt the size of the nodes to the number of genes
#  if(GeneNumber < 100)
#  {
#    Size = 11 - round(GeneNumber * 10 / 100)
#  }
  
#  if(GeneNumber >= 100)
#  {
#    Size = 1
#  }
 
#  for(i in 1:length(NodeLabel))
#  {
    ## If the node is a parent gene
#    if(sum(parentGeneList == NodeLabel[i]) > 0)
#      {
#        NodeSize = c(NodeSize, 2 * Size)	
 #       NodeColor = c(NodeColor, parentColor)
        #LabelPos = c(LabelPos, (2 * Size)/15)
 #       if(parentgeneNames == FALSE)
 #         {WritenNodeLabel[i] = ""}			
 #     }else{       
        ## If the node is a target gene
 #       if(sum(targetGeneList == NodeLabel[i]) > 0)
 #         {
 #           NodeSize = c(NodeSize, Size)
 #           NodeColor  = c(NodeColor, targetColor)
            #LabelPos = c(LabelPos, (2 * Size)/15)
 #           if(targetgeneNames == FALSE)
 #             {WritenNodeLabel[i] = ""}			
 #         }
 #     }
 # }

########################################################################
#Modified: 02/24/2012
#Last modification: 02/24/2012
#By SERVILLO David
########################################################################
#***************************************************************************************************
 
  for(i in 1:length(NodeLabel))
  {
    ## If the node is a parent gene
    if(sum(parentGeneList == NodeLabel[i]) > 0)
      {
	NodeSize = c(NodeSize, SizeParent)
        NodeColor = c(NodeColor, parentColor)
	        if(parentgeneNames == FALSE)
          {WritenNodeLabel[i] = ""}			
      }else{
        ## If the node is a target gene
        if(sum(targetGeneList == NodeLabel[i]) > 0)
          {
	    NodeSize = c(NodeSize, SizeTarget)
            NodeColor = c(NodeColor, targetColor)
            if(targetgeneNames == FALSE)
              {WritenNodeLabel[i] = ""}			

          }
      }
  }
#***************************************************************************************************

  ## One subnetwork is traced for each different CPstart in the data
  CPstartList = sort(unique(ARTIVAnet$CPstart))

  if(onepage){
    par(mfrow = c(1,length(unique(ARTIVAnet$CPstart))))
  }#else{
    #par(mfrow = c(1,1))
  #}	    
  # Representation of the global network, with all selected edges (according to T value)
#  MainTitle = paste("GLOBAL regulatory network\n", "Edge probability threshold =", Threshold)
#  PlotFunction(SubNetwork, NodeSize, NodeColor, NodeLabel, LabelPos, NodeCoord, MainTitle)
  
  # Representations of sub-networks according to temporal phases.
  for(i in 1:length(CPstartList))
    {
      CurrentCPstart = CPstartList[i]
      
      ## Delete Edges which have CPstart > CurrentCPstart or CPend < CPstart
      EdgesToDelete = unique(c(which(ARTIVAnet$CPstart > CurrentCPstart), which(ARTIVAnet$CPend < CurrentCPstart), which(ARTIVAnet$PostProb < edgesThreshold)))
      SubNetwork = delete.edges(GlobalNetwork, E(GlobalNetwork)[EdgesToDelete - 1])
  
      if(i == length(CPstartList))
  	{
          MainTitle = paste("Sub-network #", i,"\n( time point", CurrentCPstart, "to ", max(ARTIVAnet$CPend), ")")
  	}
      else
  	{
          MainTitle = paste("Sub-network #", i,"\n( time point", CurrentCPstart, "to ", CPstartList[i+1] - 1, ")") 
  	}

      PlotFunction(SubNetwork, NodeSize, NodeColor, WritenNodeLabel, LabelPos = 0, NodeCoord, MainTitle)  
      if(!onepage | i==1){
        legend(par("usr")[1],par("usr")[3],c("Positive interaction", "Negative interaction"),lty= c(1,2))
      } 	
      ## End of for() CPstartList	
    }
  
# End of function TraceNetworks()
}
