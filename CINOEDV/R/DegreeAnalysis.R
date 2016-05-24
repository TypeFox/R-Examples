DegreeAnalysis <-
function(Vertices,Edges,SaveFileName=""){
  # Degree Analysis of real vertices
  #
  # input
  #    Vertices: vertices for network construction
  #    Edges: Edges for network construction
  #    SaveFileName: basic file name for saving figure.
  #
  # output
  #    Degrees: Degrees of real vertices
  #
  # Junliang Shang
  # 3.30/2014
  
  #library(igraph)
  
  gdf <- graph.data.frame(Edges,directed=F,Vertices)
  gdfDegree <- degree(gdf)
  gdfDegree <- gdfDegree[1:length(Vertices[Vertices[,3]=="1",1])]

  gdfDegreeName <- names(gdfDegree)
  TgdfDegree <- rep(0,length(gdfDegreeName))
  for (i in 1:length(gdfDegreeName)){
    TgdfDegree[i] <- gdfDegree[[i]]
  }
  Degrees <- cbind(gdfDegreeName,TgdfDegree)
  Degrees <- Degrees[order(as.numeric(Degrees[,2]),decreasing=T),]
  colnames(Degrees) <- c("SNP","Degree")
  
  dev.new()
  set.seed(8)
  hist(gdfDegree,col="gray",border="yellow",xlab="Degree",main=NA)
  title(main=list("Histogram of real vertice degree",font=3))
  
  # savePlot(filename=paste(SaveFileName,"_DA",sep=""),type="tiff")
  
  list(Degrees=Degrees)
}
