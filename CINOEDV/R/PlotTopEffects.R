PlotTopEffects <-
function(Vertices,Top=20,SaveFileName=""){
  # plot top n effects with their corresponding SNPs or SNP-combinations.
  # Effects can be respectively considered as independent effects and combination
  # effects.
  # Independent Effect: Effect that only the SNP or SNP-combination has.
  # Combination Effect: Effect is the addition of all effects of SNP-combinaton
  # and its sub-combinations.
  #
  # input
  #    Vertices: Vertices for network construction
  #    Top: top n effects
  #    SaveFileName: basic file name for saving figure.
  #
  # output
  #    TopEffect: Independent Effects of SNPs or SNP-combinations
  #    CombinationEffect: Combination Effects of SNPs or SNP-combinations
  #
  # Junliang Shang
  # 3.28/2014
  
  if (!is.numeric(Top)){
    Top <- 20
  }
  
  if (nrow(Vertices)<Top){
    Top <- nrow(Vertices)
  }
  
  dev.new()
  par(mar=c(5,14,4,2))
  set.seed(8)
  layout(matrix(1:2, 1,2))
  
  # Independent Effect
  Vertices <- Vertices[order(as.numeric(Vertices[,2]),decreasing=T),]
  dim(Vertices) <- c(length(Vertices)/3,3)
  Values <- as.numeric(Vertices[Top:1,2])
  names(Values) <- Vertices[Top:1,1]
  barplot(Values,col=rainbow(Top),las=1,horiz=TRUE,xlab="Independent Effect")
  title(main=list(paste("Top ",Top," Effects",sep=""),font=3))
  
  # Combination Effect
  NumVertices <- nrow(Vertices)
  Tag <- rep(0,NumVertices)
  Vertices <- cbind(Vertices,Tag)
  
  for (i in 1:NumVertices){
    
    if (i==1){
      Vertices[i,4] <- "1"
      for (j in (i+1):NumVertices){
        if ((as.numeric(Vertices[i,3])>as.numeric(Vertices[j,3]))&&
              (length(strsplit(paste(Vertices[i,1]," "),Vertices[j,1])[[1]])>1)){
          Vertices[i,2] <- as.character(as.numeric(Vertices[i,2])+
                                          as.numeric(Vertices[j,2]))
        }
      } 
    }
    
    if (i==NumVertices){
      Vertices[i,4] <- "1"
      for (j in 1:(i-1)){         
        if ((as.numeric(Vertices[j,3])>as.numeric(Vertices[i,3]))&&
              (length(strsplit(paste(Vertices[j,1]," "),Vertices[i,1])[[1]])>1)){
          Vertices[i,4] <- "0"
          break
        }
        if ((as.numeric(Vertices[i,3])>as.numeric(Vertices[j,3]))&&
              (length(strsplit(paste(Vertices[i,1]," "),Vertices[j,1])[[1]])>1)){
          Vertices[i,2] <- as.character(as.numeric(Vertices[i,2])+
                                          as.numeric(Vertices[j,2]))
          Vertices[j,4] <- "0"
        }
      }
    }
    
    if ((i>1)&&(i<NumVertices)){
      Vertices[i,4] <- "1"
      for (j in (i+1):NumVertices){
        if ((as.numeric(Vertices[i,3])>as.numeric(Vertices[j,3]))&&
              (length(strsplit(paste(Vertices[i,1]," "),Vertices[j,1])[[1]])>1)){
          Vertices[i,2] <- as.character(as.numeric(Vertices[i,2])+
                                          as.numeric(Vertices[j,2]))
        }
      }
      for (j in 1:(i-1)){         
        if ((as.numeric(Vertices[j,3])>as.numeric(Vertices[i,3]))&&
              (length(strsplit(paste(Vertices[j,1]," "),Vertices[i,1])[[1]])>1)){
          Vertices[i,4] <- "0"
          break
        }
        if ((as.numeric(Vertices[i,3])>as.numeric(Vertices[j,3]))&&
              (length(strsplit(paste(Vertices[i,1]," "),Vertices[j,1])[[1]])>1)){
          Vertices[i,2] <- as.character(as.numeric(Vertices[i,2])+
                                          as.numeric(Vertices[j,2]))
          Vertices[j,4] <- "0"
        }
      }
    }
    
    if (sum(as.numeric(Vertices[,4]))==Top){
      break
    }
  }
  
  Vertices <- Vertices[Vertices[,4]=="1",1:3]
  
  Top <- nrow(Vertices)
  Vertices <- Vertices[order(as.numeric(Vertices[,2]),decreasing=T),]
  Values2 <- as.numeric(Vertices[Top:1,2])
  names(Values2) <- Vertices[Top:1,1]
  barplot(Values2,col=rainbow(Top),las=1,horiz=TRUE,xlab="Combination Effect")
  title(main=list(paste("Top ",Top," Effects",sep=""),font=3))
  
  # save and return
  # savePlot(filename=paste(SaveFileName,"_TE",sep=""),type="tiff")
  
  list(TopEffect=Values,CombinationEffect=Values2)  
}
