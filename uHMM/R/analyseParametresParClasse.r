#' @title Generate graphics by states
#' @description This function is used by the \code{\link{uHMMinterface}} to create plots and boxplots by state, and state sequencing after classification, modeling and prediction steps.
#' @param data matrix or dataframe.
#' @param stateSeq state sequencing obtained by classification, modeling or prediction step.
#' @param step character, one of the following : "classification", "modeling" or "prediction".
#' @param period numeric vector of observation times.
#' @param directory path of the directory in which a sub-directory will be created to save graphics.
#' @param uHMMinterface logical : whether the function is used via the uHMM interface.
#' @param graphicFrame frame in which the state sequencing should be displayed.
#' @param hscaleGraphicFrame the hscale parameter value of the tkplot function used to create the graphic frame.
#' @param vscaleGraphicFrame the vscale parameter value of the tkplot function used to create the graphic frame.
#' @param tm a one row dataframe containing text to display in the interface.
#' @return Plots and boxplots for each variable of the dataset and the state sequencing are saved in the diagram subdirectory.
#' @return Summary by cluster is saved in the table subdirectory.
#' @importFrom graphics axis par

.graphicsByState<-function(data,stateSeq,step=NULL,period,directory,uHMMinterface=FALSE,hscaleGraphicFrame=1.2,vscaleGraphicFrame=1.2,tm,graphicFrame){
  
  
  # Fonction pour mettre premiere lettre d'une chaine en majuscule
  .simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")
  }
  
 outFigures<-paste(directory,tm$diagramsRepertory,"/",sep="")
 outTables<-paste(directory,tm$tablesRepertory,"/",sep="")
  
  #Initialisation des constantes et variables
  data2=data[,-c(which(colnames(data)=="Dates"), which(colnames(data)=="Hours"))]
  nbVar=dim(data2)[2]
  nbCluster=max(stateSeq,na.rm=TRUE)
  coul=2:(nbCluster+1)
  
  #print(paste("gap",uHMMenv$gap))

  tableLine=1
  summaryTable=matrix("",2*(nbCluster-1)*dim(data2)[2],7)
  
  #Boxplot et representation de chque parametre par cluster
  for(a in 1:nbVar){    #Pour chaque parametre
    x=data2[,a];
    for (i in 2:nbCluster){ #Pour chaque Cluster dans ce parametre
      
      #Completion du tableau des statistiques de base
      summaryTable[tableLine,1]=colnames(data2)[a]
      summaryTable[tableLine,2]=paste("Cluster_",i-1,sep="")
      if(length(summary(x[stateSeq==i]))==7){
        summaryTable[tableLine+1,]=summary(x[stateSeq==i])
      }
      else{
        summaryTable[tableLine+1,]=c(summary(x[stateSeq==i]),0)
      }

      #Creation d'une liste pour la creation des graphes
      if (i==2){
        listePar=list(x[stateSeq==i])
      }
      else{
        listePar=c(listePar,list(x[stateSeq==i]))
      }   
      tableLine=tableLine+2
      
    } 

    #Boite a moustache de chaque parametre pour chaque cluster
    figure=paste(outFigures,paste("boxplot",colnames(data2)[a],".jpg",sep=""),sep="")
    jpeg(figure,700,700)
    boxplot(listePar,col=coul,main=colnames(data2)[a],xaxt = "n",cex.axis=1.5,cex.lab=1.5,cex=1.5,xlab="Cluster",ylab="Level")
    axis(1, at=1:(nbCluster-1), labels=c(1:(nbCluster-1)), tick = TRUE, col = "black",cex.axis=1.5)
    dev.off()
    
    #Representation de chaque parametre avec la projection des clusters
    figure=paste(outFigures,paste("Representation_",colnames(data2)[a],".jpg",sep=""),sep="")
    jpeg(figure,700,700)
    plot(data2[,a]~chron(period),col=stateSeq,main=colnames(data2)[a],xlab="Time",ylab="Level",cex.axis=1.5,cex.lab=1.5,cex=1.5)
    dev.off()
  }

  seqTitle<-NULL
  if (step=="classification"){seqTitle<-tm$seqClassifTitle}
  if (step=="modeling"){seqTitle<-tm$seqModelTitle}
  if (step=="prediction"){seqTitle<-tm$seqPredictionTitle}  
  
  #Representation du sequencement des classes pour la periode choisie
  figure=paste(outFigures,tm$seqClassifFileTitle,sep="")
  jpeg(figure,700,700)
  plot(stateSeq~chron(period),col=stateSeq,main=seqTitle,xlab="Time",ylab="Cluster",cex.axis=1.5,cex.lab=1.5,cex=1.5,yaxt="n")
  axis(2, at=1:nbCluster, labels=c("NA",1:(nbCluster-1)), tick = TRUE, col = "black",cex.axis=1.5)
  dev.off()

  if (uHMMinterface==TRUE){
    # Affichage du sequencement des classes pour la periode choisie
    seqClassPlot<-tkrplot(graphicFrame,hscale=hscaleGraphicFrame,vscale=vscaleGraphicFrame,function(){
      par(bg = "white")
      plot(stateSeq~chron(period),col=stateSeq,main=seqTitle,xlab="Time",ylab="Cluster",cex.axis=1,cex.lab=1,cex=0.8,yaxt="n")
      axis(2, at=1:nbCluster, labels=c("NA",1:(nbCluster-1)), tick = TRUE, col = "black",cex.axis=1)
    })
    tkgrid(seqClassPlot,row=0,column=1,sticky="w")
  }
  

  #Modification min max du boxplot
  for(a in 1:nbVar){ #Pour chaque parametre
    x=data2[,a];
    qu3<-c() # 3e Quartile 
    qu1<-c() # 1er Quartile
    for (i in 2:nbCluster){ #Pour chaque Cluster dans ce parametre
      
      if (i==2){
        listePar=list(x[stateSeq==i])
      }
      else{
        listePar=c(listePar,list(x[stateSeq==i]))
      }
      
      qu3<-c(qu3,summary(x[stateSeq==i])[5]) 
      qu1<-c(qu1,summary(x[stateSeq==i])[2])  
    } 
  
    #Boite a moustache de chaque parametre pour chaque classe
    figure=paste(outFigures,paste("boxplot",colnames(data2)[a],"_Reduit.jpg",sep=""),sep="")
    jpeg(figure,700,700)
    boxplot(listePar,col=coul,main=colnames(data2)[a],xaxt = "n",ylim=c((min(qu1,na.rm=TRUE)-min(qu1,na.rm=TRUE)*0.5),max(qu3,na.rm=TRUE)+(max(qu3,na.rm=TRUE)*1.5)),cex.axis=1.5,cex.lab=1.5,cex=1.5,xlab="Cluster",ylab="Level");
    axis(1, at=1:(nbCluster-1), labels=c(1:(nbCluster-1)), tick = TRUE, col = "black",cex.axis=1.5)
    dev.off()
  }

  colnames(summaryTable)=c("Min","Q1","Median","Mean","Q3","Max","NA")
  #Sauvegarde
  write.table(summaryTable,file=paste(outTables,"summaryTable",.simpleCap(step),".xls",sep=""),row.names=FALSE)
  
}

    