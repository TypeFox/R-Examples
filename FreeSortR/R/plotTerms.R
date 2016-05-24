
#####################################################################
#_____________    function plotTerms    ______________________
#                 function plotTerms
#####################################################################
plotTerms<-function(MatTerms,ResMds,dim=c(1,2),type="correl",add=TRUE){
  
  
  if (!class(ResMds)=="SortingMds"){
    return("ResMds is not an object of class SortingMds.")
  } else {
    Config<-ResMds@Config
    if (nrow(MatTerms)!=nrow(Config)){
      return("Check the dimension of MatTerms array.")
    } else {
      dim1<-dim[1]
      dim2<-dim[2]
      xtex=paste("Dim ",dim1," (",round(100*ResMds@Percent[dim1],2)," %)",sep="")
      ytex=paste("Dim ",dim2," (",round(100*ResMds@Percent[dim2],2)," %)",sep="")
      
      if (type=="correl"){
        
        # correlations between terms and dimensions of Mds
        Correl<-cor(MatTerms,Config)
        
        
        xtex=paste("Dim ",dim1,sep="")
        ytex=paste("Dim ",dim2,sep="")
        
        plot(Correl[,dim1],Correl[,dim2],type="n",asp=1,xlim=c(-1,1),ylim=c(-1,1),xlab=xtex,ylab=ytex,lwd=5,pch=16)
        text(Correl[,dim1],Correl[,dim2],colnames(MatTerms))
        abline(h=0,v=0)
        symbols(0,0,circles=1,inches=FALSE,add=TRUE)
        arrows(0,0,Correl[,dim1],Correl[,dim2],length=0.10)
        
        rownames(Correl)<-colnames(MatTerms)
        
        return(Correl)
        
      } else {
        
        #barycentric representation
        
        weights<-sweep(MatTerms,MARGIN=2,STATS=apply(MatTerms,2,sum),FUN="/")
        Baryc<-t(weights)%*%Config
        
        xtex=paste("Dim ",dim1," (",round(100*ResMds@Percent[dim1],2)," %)",sep="")
        ytex=paste("Dim ",dim2," (",round(100*ResMds@Percent[dim2],2)," %)",sep="")
        
        
        if (add==TRUE){
          
          xlim<-c(min(min(Baryc[,dim1]),min(Config[,dim1])),max(max(Baryc[,dim1]),max(Config[,dim1])))
          ylim<-c(min(min(Baryc[,dim2]),min(Config[,dim2])),max(max(Baryc[,dim2]),max(Config[,dim2])))
          plot(Baryc[,dim1],Baryc[,dim2],type="n",xlab=xtex,ylab=ytex,xlim=xlim,ylim=ylim)
          text(Baryc[,dim1],Baryc[,dim2],colnames(MatTerms))
          abline(h=0,v=0)
          text(Config[,dim[1]],Config[,dim[2]],ResMds@LabStim,col="blue")
          
        } else {
          
          plot(Baryc[,dim1],Baryc[,dim2],type="n",xlab=xtex,ylab=ytex)
          abline(h=0,v=0)
          text(Config[,dim[1]],Config[,dim[2]],ResMds@LabStim,col="blue")
        }
        
        
        
        rownames(Baryc)<-colnames(MatTerms)
        
        return(Baryc)
        
      }
      
    }
  }
}
