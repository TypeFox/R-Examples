drawSpec <-function(X, startP=-1, endP=-1, groupLabel=NULL, useLog=-1,
     highBound=-1, lowBound=-1, xlab=NULL, ylab=NULL, main=NULL,
     nAxisPos=4, offside=0)
{
    groupLabel_name=groupLabel;
    X=as.data.frame(X);
    colnames(X)=c(1:ncol(X));
    X=as.matrix(X);
    if (highBound!=-1){
        for (i in 1:nrow(X)){
            myIndex=which(X[i,]>highBound);
            X[i,myIndex]=highBound;
        }        
    }
    
    if (lowBound!=-1){
        for (i in 1:nrow(X)){
            myIndex=which(X[i,]<lowBound);
            X[i,myIndex]=lowBound;
        }        
    }
        
    if (is.null(groupLabel)){
        groupLabel=c(1:nrow(X));
        groupLabel=as.factor(groupLabel);
    }else{
        levels(groupLabel) =c(1:length(levels(groupLabel)))
    }
    
    if (startP==-1) startP=1;
    if (endP==-1) endP=ncol(X);
    
    if (is.null(xlab)){xlab='index' }
    if (is.null(ylab)){ylab='intensity' }
    if (is.null(main)){main=paste(' ',startP+offside,'-',endP+offside)}
    
    GraphRange<-c(startP:endP); 
    
    yn<-X[,GraphRange];
    if (useLog!=-1) yn=log(yn);
        
    plot(yn[1,],ylim=c(min(yn),max(yn)),type="n",
    ylab=ylab,
    xlab=xlab,
    main=main,
    xaxt="n" 
    )
    tempVal =trunc(length(GraphRange)/nAxisPos);
    xPos=c(0:nAxisPos) * tempVal; 
    axis(1,at=xPos,labels=xPos+startP+offside);
        
    for(i in 1:length(levels(groupLabel))){
      groupLabelIdx=which(groupLabel==levels(groupLabel)[i]);
      for (j in 1:length(groupLabelIdx)){
        lines(yn[groupLabelIdx[j],],col=as.integer(levels(groupLabel)[i]))        
      }
    }

     if (!is.null(groupLabel_name)){
         legendPos="topleft";
         legend(legendPos,levels(groupLabel_name),col=as.integer(levels(groupLabel)),
         text.col="black",pch=c(19,19),bg='gray90')
     }
}
