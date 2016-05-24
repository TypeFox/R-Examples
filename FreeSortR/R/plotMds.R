
#####################################################################
#_________________    function plotMds          ______________________
#               function plotMds for class SortingMds
#####################################################################
plotMds <- function(ResMds,dim=c(1,2),ellipse=FALSE,proba=0.90,col=NULL){ 
  
  if (!class(ResMds)=="SortingMds"){
    return("This is not an object of class SortingMds")
  }  
  else
  {  
    if (ellipse==FALSE){    
      dim1<-dim[1]
      dim2<-dim[2]
      xtex=paste("Dim ",dim1," (",round(100*ResMds@Percent[dim1],2)," %)",sep="")
      ytex=paste("Dim ",dim2," (",round(100*ResMds@Percent[dim2],2)," %)",sep="")
      
      plot(ResMds@Config[,dim1],ResMds@Config[,dim2],type="n",xlab=xtex,ylab=ytex)
      text(ResMds@Config[,dim1],ResMds@Config[,dim2],labels=ResMds@LabStim,col=col)
    }else{
      if (length(ResMds@ResBoot)!=ResMds@nstimuli){
        cat("Problem of dimension of bootstrap results.")
      } else {
        ResBoot2<-ResMds@ResBoot
        nstim<-length(ResBoot2)
        ndim<-dim(ResBoot2[[1]])[2]
        nboot<-dim(ResBoot2[[1]])[1]
        
        ResBoot<-vector("list")
        for (stim in 1:nstim){
          Mat<-ResBoot2[[stim]]
          ResBoot[[stim]]<-Mat[,dim]
        }
        
        #library(ellipse)
        
        CentreE=matrix(0,nstim,2)  
        for (stim in 1:nstim){
          CentreE[stim,]=colMeans(ResBoot[[stim]])
        }   
        
        bid1<-matrix(0,nstim,2)
        bid2<-matrix(0,nstim,2)
        for (stim in 1:nstim){
          bid1[stim,]<-apply(ResBoot[[stim]],MARGIN=2,FUN=max)
          bid2[stim,]<-apply(ResBoot[[stim]],MARGIN=2,FUN=min)
        }   
        
        max1<-max(bid1[,1])
        max2<-max(bid1[,2])
        min1<-min(bid2[,1])
        min2<-min(bid2[,2])
        
        dim1<-dim[1]
        dim2<-dim[2]
        xtex=paste("Dim ",dim1," (",round(100*ResMds@Percent[dim1],2)," %)",sep="")
        ytex=paste("Dim ",dim2," (",round(100*ResMds@Percent[dim2],2)," %)",sep="")
        
        
        plot(CentreE[,1],CentreE[,2],type="p",xlim=c(min1,max1),ylim=c(min2,max2),xlab=xtex,ylab=ytex) 
        text(CentreE[,1],CentreE[,2],label=ResMds@LabStim)  
        
        for (stim in 1:nstim){
          tabstim=ResBoot[[stim]]
          covstim<-cov(tabstim)
          res<-ellipse(covstim[1:2,1:2],centre=CentreE[stim,],level=proba)
          lines(res)
        }
      }    
    }
  }
}
