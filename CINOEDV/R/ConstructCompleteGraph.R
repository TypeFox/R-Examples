ConstructCompleteGraph <-
function(Vertices,Edges,BaseSize=6,SaveFileName=""){
  # Construct Complete Graph by using all vertices and edges
  #
  # input
  #    vertices:  Vertices for network construction
  #    Edges:  Edges for network construction
  #    BaseSize: Basic size of vertices
  #    SaveFileName: basic file name for saving figure.
  #
  # Junliang Shang
  # 3.28/2014
  
  #library(igraph)
  
  if((nrow(Edges)==1) && (Edges[1,1]==0) && (Edges[1,2]==0)){
    if((nrow(Vertices)==1) && (Vertices[1,1]==0) && (Vertices[1,2]==0) && (Vertices[1,3]==0)) {
      cat("  Both Vertices and Edges are empty !\n")
    }else
    {
      
      Edges <- c(Vertices[,1],rev(Vertices[,1]))
      dim(Edges) <- c(length(Edges)/2,2)
      gdf <- graph.data.frame(Edges,directed=F,Vertices)
      
      dev.new()
      set.seed(8)
      par(mar = c(0, 0, 0, 0))
      
      # Vertice Colour
      Colour <- rep("red", nrow(Vertices))
      
      # Vertice Label
      Label <- Vertices[,1]
      dim(Label)=c(length(Label),1)
      
      # Vertice Size
      if (nrow(Vertices)==1){
        Size <- BaseSize
      }else
      {
        Size <- as.numeric(Vertices[,2])
        dim(Size) <- c(length(Size),1)
        minSize <- min(Size)
        maxSize <- max(Size)
        for (i in 1:length(Size)){
          Size[i,1] <- (((Size[i,1]-minSize)/(maxSize-minSize))*2+1)*BaseSize
        }
      }
      
      plot(gdf,layout=layout.circle, vertex.label=Label,vertex.color=Colour,vertex.size=Size,
           vertex.frame.color=Colour,edge.lty=0)
    }
    
  }else
  {
    
    gdf <- graph.data.frame(Edges,directed=F,Vertices)
    dev.new()
    set.seed(8)
    par(mar = c(0, 0, 0, 0))
    
    # Vertice Colour
    Colour <- rep("red", nrow(Vertices))
    Colour[Vertices[,3]=="1"] <- "red"
    Colour[Vertices[,3]=="2"] <- "gray"
    Colour[Vertices[,3]=="3"] <- "darkslategray1"
    Colour[Vertices[,3]=="4"] <- "darkkhaki"
    Colour[Vertices[,3]=="5"] <- "darksalmon"
    
    # Vertice Label
    Label <- Vertices[,1]
    dim(Label)=c(length(Label),1)
    Label[Vertices[,3]!="1",] <- NA
    
    # Vertice Size
    Size <- as.numeric(Vertices[,2])
    dim(Size) <- c(length(Size),1)
    minSize <- min(Size)
    maxSize <- max(Size)
    for (i in 1:length(Size)){
      Size[i,1] <- (((Size[i,1]-minSize)/(maxSize-minSize))*2+1)*BaseSize
    }
    
    # plot
    plot(gdf, layout=layout.fruchterman.reingold, vertex.label=Label, 
         vertex.color=Colour,vertex.size=Size,vertex.frame.color=Colour)    
  }
  
  # savePlot(filename=paste(SaveFileName,"_CG",sep=""),type="tiff")
}
