FNPre<-function(file,indicator=c("RA","AA","Jaccard"),threshold=0.1, output="FalseNegativePreResult-ppiPre.csv") 
{
  edges<-read.csv(file=file,header=TRUE,sep=",")
	graph<-graph.data.frame(edges,directed=FALSE)
	nodes<-get.vertex.attribute(graph,"name")
	cache_jaccard<-similarity.jaccard(graph) #jaccard similarity matrix
	x<-array(1:(length(nodes)*(length(nodes)-1)/2))	
	TopSims<-data.frame(protein1="",protein2="",Jaccard=x,AA=x,RA=x,label=0) 
	TopSims$protein1<-as.character(TopSims$protein1) 
	TopSims$protein2<-as.character(TopSims$protein2) 
	i<-1
	y<-1
	
	for(n in 1:length(nodes))
	{
		if(n!=length(nodes))
		{	for(m in (n+1):length(nodes))
			{
				TopSims[[1]][i]<-nodes[n]
				TopSims[[2]][i]<-nodes[m]
				i<-i+1
			}
		}
	}
        
       for(n in 1:length(x)) 
	{
		for(m in 1:length(edges[[1]]))
		{
			if(as.character(TopSims[[1]][n])==as.character(edges[[1]][m]))
			{
				if(as.character(TopSims[[2]][n])==as.character(edges[[2]][m]))
				{
					TopSims[[6]][n]<-1
					break
				}
			}
			else if(as.character(TopSims[[1]][n])==as.character(edges[[2]][m]))
			{
				if(as.character(TopSims[[2]][n])==as.character(edges[[1]][m]))
				{
					TopSims[[6]][n]<-1
					break
				}
			}
		}
	}

     	ifJaccard<-match("Jaccard",indicator)   #if Jaccard similarity is used
	ifAA<-match("AA",indicator)
	ifRA<-match("RA",indicator)
	
	if(is.na(ifJaccard))
		ifJaccard<-"NA"
	if(is.na(ifAA))
		ifAA<-"NA"
	if(is.na(ifRA))
		ifRA<-"NA"
	

        if(ifJaccard!="NA")
	{        
		i<-1
		for(n in 1:length(nodes)) 
		{
			if(n!=length(nodes))
			{	for(m in (n+1):length(nodes))
				{
					TopSims[[3]][i]<-cache_jaccard[n,m]
					i<-i+1
				}
			}
		}
	}

	if(ifAA!="NA")
	{
		for(i in 1:length(x)) # calculate AA similarity between two nodes
		{
			TopSims[[4]][i]<-AASim(as.character(TopSims[[1]][i]),as.character(TopSims[[2]][i]),graph)
		}
	}

	if(ifRA!="NA")
	{
		for(i in 1:length(x)) # calculate RA similarity between two nodes
		{
			TopSims[[5]][i]<-RASim(as.character(TopSims[[1]][i]),as.character(TopSims[[2]][i]),graph)
		}
	}

	#build negative data set based on "indicator" and "threshold"
	threshold_number=length(x)*threshold
	simi_result<-data.frame(protein1="",protein2="",Jaccard=x,AA=x,RA=x,label=0) 
	simi_result$protein1<-as.character(simi_result$protein1)
	simi_result$protein2<-as.character(simi_result$protein2) 
	
	if(ifJaccard!="NA"&&ifAA=="NA"&&ifRA=="NA")  #indicator="Jaccard"
	{
		simi_jaccard<-TopSims[order(TopSims[,3],decreasing=TRUE),]  #sort rows in TopSims based on their Jaccard value
		simi_jaccard_neg <- simi_jaccard[simi_jaccard[[6]]==0,] #only use negative data
		for(n in 1:threshold_number)
		{
			simi_result[[1]][n]<-simi_jaccard_neg [[1]][n]
			simi_result[[2]][n]<-simi_jaccard_neg [[2]][n]
			simi_result[[3]][n]<-simi_jaccard_neg [[3]][n]
			simi_result[[4]][n]<-simi_jaccard_neg [[4]][n]
			simi_result[[5]][n]<-simi_jaccard_neg [[5]][n]
			simi_result[[6]][n]<-simi_jaccard_neg [[6]][n]
		}
		y<-threshold_number
		simi_result$AA<-0
		simi_result$RA<-0
	}

	if(ifJaccard=="NA"&&ifAA!="NA"&&ifRA=="NA")  #indicator="AA"
	{
		simi_AA<-TopSims[order(TopSims[,4],decreasing=TRUE),]  #sort rows in TopSims based on their AA value
		simi_AA_neg <- simi_AA[simi_AA[[6]]==0,] #only use negative data
		for(n in 1:threshold_number)
		{
			simi_result[[1]][n]<-simi_AA_neg [[1]][n]
			simi_result[[2]][n]<-simi_AA_neg [[2]][n]
			simi_result[[3]][n]<-simi_AA_neg [[3]][n]
			simi_result[[4]][n]<-simi_AA_neg [[4]][n]
			simi_result[[5]][n]<-simi_AA_neg [[5]][n]
			simi_result[[6]][n]<-simi_AA_neg [[6]][n]
		}
		y<-threshold_number
		simi_result$Jaccard<-0
		simi_result$RA<-0
	}

	if(ifJaccard=="NA"&&ifAA=="NA"&&ifRA!="NA")  #indicator="RA"
	{
		simi_RA<-TopSims[order(TopSims[,5],decreasing=TRUE),]  #sort rows in TopSims based on their RA value
		simi_RA_neg <- simi_RA[simi_RA[[6]]==0,] #only use negative data
		for(n in 1:threshold_number) #select negative data, the number of data entries= number of positive data entries * threshold
		{
			simi_result[[1]][n]<-simi_RA_neg [[1]][n]
			simi_result[[2]][n]<-simi_RA_neg [[2]][n]
			simi_result[[3]][n]<-simi_RA_neg [[3]][n]
			simi_result[[4]][n]<-simi_RA_neg [[4]][n]
			simi_result[[5]][n]<-simi_RA_neg [[5]][n]
			simi_result[[6]][n]<-simi_RA_neg [[6]][n]
		}
		y<-threshold_number
		simi_result$Jaccard<-0
		simi_result$AA<-0
	}

	if(ifJaccard!="NA"&&ifAA!="NA"&&ifRA=="NA") #indicator=c("Jaccard","AA")
	{
		simi_jaccard<-TopSims[order(TopSims[,3],decreasing=TRUE),]
		simi_jaccard_neg <- simi_jaccard[simi_jaccard[[6]]==0,] #only use negative data
		simi_AA<-TopSims[order(TopSims[,4],decreasing=TRUE),]
		simi_AA_neg <- simi_AA[simi_AA[[6]]==0,] #only use negative data
		y<-1
		for(n in 1:threshold_number)
    		{
			for(m in 1:threshold_number)
   			{
				if(simi_jaccard_neg[[1]][n]==simi_AA_neg[[1]][m]&&simi_jaccard_neg[[2]][n]==simi_AA_neg[[2]][m])
				{
					simi_result[[1]][y]<-simi_jaccard_neg[[1]][n]
					simi_result[[2]][y]<-simi_jaccard_neg[[2]][n]
					simi_result[[3]][y]<-simi_jaccard_neg[[3]][n]
					simi_result[[4]][y]<-simi_jaccard_neg[[4]][n]
					simi_result[[5]][y]<-simi_jaccard_neg[[5]][n]
					simi_result[[6]][y]<-simi_jaccard_neg[[6]][n]
					y <- y+1
					break
				}
			}
		}
		y <- y-1
		simi_result$RA <- 0
	}

	if(ifJaccard!="NA"&&ifAA=="NA"&&ifRA!="NA")   #indicator=c("Jaccard",RA)
	{
		simi_jaccard<-TopSims[order(TopSims[,3],decreasing=TRUE),]
		simi_jaccard_neg <- simi_jaccard[simi_jaccard[[6]]==0,] #only use negative data

		simi_RA<-TopSims[order(TopSims[,5],decreasing=TRUE),]
		simi_RA_neg <- simi_RA[simi_RA[[6]]==0,] #only use negative data
		y<-1
		for(n in 1:threshold_number)
    		{
			for(m in 1:threshold_number)
   			{
				if(simi_jaccard_neg[[1]][n]==simi_RA_neg[[1]][m]&&simi_jaccard_neg[[2]][n]==simi_RA_neg[[2]][m])
				{
					simi_result[[1]][y]<-simi_jaccard_neg[[1]][n]
					simi_result[[2]][y]<-simi_jaccard_neg[[2]][n]
					simi_result[[3]][y]<-simi_jaccard_neg[[3]][n]
					simi_result[[4]][y]<-simi_jaccard_neg[[4]][n]
					simi_result[[5]][y]<-simi_jaccard_neg[[5]][n]
					simi_result[[6]][y]<-simi_jaccard_neg[[6]][n]
					y<-y+1
					break
				}
			}
		}
		y<-y-1
		simi_result$AA<-0
	}
        
	if(ifJaccard=="NA"&&ifAA!="NA"&&ifRA!="NA")   #indicator=c("AA",RA)
	{
		simi_AA<-TopSims[order(TopSims[,4],decreasing=TRUE),]
		simi_AA_neg <- simi_AA[simi_AA[[6]]==0,] #only use negative data

		simi_RA<-TopSims[order(TopSims[,5],decreasing=TRUE),]
		simi_RA_neg <- simi_RA[simi_RA[[6]]==0,] #only use negative data

		y<-1
		for(n in 1:threshold_number)
    		{
			for(m in 1:threshold_number)
   			{
				if(simi_AA_neg[[1]][n]==simi_RA_neg[[1]][m]&&simi_AA_neg[[2]][n]==simi_RA_neg[[2]][m])
				{
					simi_result[[1]][y]<-simi_AA_neg[[1]][n]
					simi_result[[2]][y]<-simi_AA_neg[[2]][n]
					simi_result[[3]][y]<-simi_AA_neg[[3]][n]
					simi_result[[4]][y]<-simi_AA_neg[[4]][n]
					simi_result[[5]][y]<-simi_AA_neg[[5]][n]
					simi_result[[6]][y]<-simi_AA_neg[[6]][n]
					y<-y+1
					break
				}
			}
		}
		y<-y-1
		simi_result$Jaccard<-0
	}

	if(ifJaccard!="NA"&&ifAA!="NA"&&ifRA!="NA") #indicator=c("Jaccard","AA","RA")
	{
		simi_jaccard<-TopSims[order(TopSims[,3],decreasing=TRUE),]
		simi_jaccard_neg <- simi_jaccard[simi_jaccard[[6]]==0,] #only use negative data

		simi_AA<-TopSims[order(TopSims[,4],decreasing=TRUE),]
		simi_AA_neg <- simi_AA[simi_AA[[6]]==0,] #only use negative data

		simi_RA<-TopSims[order(TopSims[,5],decreasing=TRUE),]
		simi_RA_neg <- simi_RA[simi_RA[[6]]==0,] #only use negative data

		y<-1
		for(n in 1:threshold_number)
    		{
			for(m in 1:threshold_number)
   			{
				if(simi_jaccard_neg[[1]][n]==simi_AA_neg[[1]][m]&&simi_jaccard_neg[[2]][n]==simi_AA_neg[[2]][m])
				{
					
					for(p in 1:threshold_number)
					{
						if(simi_jaccard_neg[[1]][n]==simi_RA_neg[[1]][p]&&simi_jaccard_neg[[2]][n]==simi_RA_neg[[2]][p])
						{
							simi_result[[1]][y]<-simi_jaccard_neg[[1]][n]
							simi_result[[2]][y]<-simi_jaccard_neg[[2]][n]
							simi_result[[3]][y]<-simi_jaccard_neg[[3]][n]
							simi_result[[4]][y]<-simi_jaccard_neg[[4]][n]
							simi_result[[5]][y]<-simi_jaccard_neg[[5]][n]
							simi_result[[6]][y]<-simi_jaccard_neg[[6]][n]
							y<-y+1
							break
						}
					}

				}
			}
		}
		y<y-1
	}
	
	z<-array(1:y)
	simi_final<-data.frame(protein1="",protein2="",Jaccard=z,AA=z,RA=z,label=0)
	simi_final$protein1<-as.character(simi_final$protein1)
	simi_final$protein2<-as.character(simi_final$protein2)
	for(i in 1:y)
	{
		simi_final[[1]][i]<-simi_result[[1]][i]
		simi_final[[2]][i]<-simi_result[[2]][i]
		simi_final[[3]][i]<-simi_result[[3]][i]
		simi_final[[4]][i]<-simi_result[[4]][i]
		simi_final[[5]][i]<-simi_result[[5]][i]
		simi_final[[6]][i]<-simi_result[[6]][i]
	}
        	
	write.csv(simi_final,file=output,row.names=FALSE) 
}



JaccardSim <- function(node1, node2,graph)    #Jaccard similarity of two nodes
{
       #if(!require("igraph")){ stop("package igraph is needed.")}
  #node1<-"1134"
  #node2<-"1147"
	nodes<-get.vertex.attribute(graph,"name")
  neighbor1 <- lapply(neighborhood(graph,1,node1), as.vector)
	neighbor2 <- lapply(neighborhood(graph,1,node2), as.vector)
	neighbor1<-neighbor1[[1]][-1] #remove the node itself
	neighbor2<-neighbor2[[1]][-1] #remove the node itself
	#neighbor1 <- neighborhood(graph,1,node1) #get neighbors of node1 including node1 itself
	#neighbor1[[1]][1]<-length(nodes)+1  #remove the node itself
	#neighbor2 <- neighborhood(graph,1,node2)
	#neighbor2[[1]][1]<-length(nodes)+2
	#insersec <- length(na.omit(match(neighbor1[[1]],neighbor2[[1]])))
	insersec <- length(na.omit(match(neighbor1,neighbor2)))
	uni <- length(neighbor1)+length(neighbor2)-insersec
	sim <- insersec/uni
	return(sim)
}

AASim<- function(node1,node2,graph) #AA similarity of two nodes
{
  #if(!require("igraph")){ stop("package igraph is needed.")}
  nodes<-get.vertex.attribute(graph,"name")	
  neighbor1 <- lapply(neighborhood(graph,1,node1), as.vector)
  neighbor2 <- lapply(neighborhood(graph,1,node2), as.vector)
  neighbor1<-neighbor1[[1]][-1] #remove the node itself
  neighbor2<-neighbor2[[1]][-1] #remove the node itself
  commons<-na.omit(match(neighbor1,neighbor2)) #common neighbours of two nodes
  sim<-0
  if(length(commons)!=0)
    for(n in 1:length(commons))
      sim<-sim+1/log10(degree(graph,v=nodes[neighbor2[commons[n]]]))
  
  return(sim)
}

RASim<-function(node1,node2,graph) #RA similarity of two nodes
{
  #if(!require("igraph")){ stop("package igraph is needed.")}   
  nodes<-get.vertex.attribute(graph,"name")	
  neighbor1 <- lapply(neighborhood(graph,1,node1), as.vector)
  neighbor2 <- lapply(neighborhood(graph,1,node2), as.vector)
  neighbor1<-neighbor1[[1]][-1] #remove the node itself
  neighbor2<-neighbor2[[1]][-1] #remove the node itself
  commons<-na.omit(match(neighbor1,neighbor2)) #common neighbours of two nodes
  sim<-0
  if(length(commons)!=0)
  {	for(n in 1:length(commons))
  {
    sim<-sim+1/degree(graph,v=nodes[neighbor2[commons[n]]])
  }
  }
  return(sim)
}

TopologicSims<-function(inputfile,outputfile="TopologicSims-ppiPre.csv", header=TRUE, sep=",") 
{
       #if(!require("igraph")){ stop("package igraph is needed.")}
	proteinname<-read.csv(file=inputfile,header=header,sep=sep)
	graph<-graph.data.frame(proteinname,directed=FALSE)
     	sims<-data.frame(protein1=proteinname[1],protein2=proteinname[2],Jaccard=0,AA=0,RA=0)
	
	#sims$protein1<-character(sims$protein1)
	#sims$protein2<-character(sims$protein2)
	
	for(i in 1:length(proteinname[[1]]))
	{
		sims[[3]][i]<-JaccardSim(as.character(sims[[1]][i]),as.character(sims[[2]][i]),graph)
		sims[[4]][i]<-AASim(as.character(sims[[1]][i]),as.character(sims[[2]][i]),graph)
		sims[[5]][i]<-RASim(as.character(sims[[1]][i]),as.character(sims[[2]][i]),graph)
		i<-i+1
	}
	write.csv(sims,file=outputfile,row.names=FALSE) 
}