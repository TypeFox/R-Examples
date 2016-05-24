graph<-function(design,data=read.table(file.choose(new=FALSE)),xlab="Measurement Times",ylab="Scores"){

  MT<-nrow(data)
  x<-1:MT
  
  if(design=="CRD"|design=="RBD"|design=="ATD"){
    plot(x,data[,2],type="n",xlab=xlab,ylab=ylab)
    points(x[data[,1]=="A"],data[,2][data[,1]=="A"],pch=1)
    points(x[data[,1]=="B"],data[,2][data[,1]=="B"],pch=16)
    a<-data[,2][data[,1]=="A"]
    b<-data[,2][data[,1]=="B"]
    MTA<-x[data[,1]=="A"]
    MTB<-x[data[,1]=="B"]
    for(it in 1:(length(a)-1)){
      lines(c(MTA[it],MTA[it+1]),c(a[it],a[it+1]),lty=2)
    }
    for(it in 1:(length(b)-1)){
      lines(c(MTB[it],MTB[it+1]),c(b[it],b[it+1]),lty=1)
    }
    legend(locator(1),lty=c(2,1),pch=c(1,16),legend=c("A","B"))
  }
  
  if(design=="AB"){
    plot(x,data[,2],xlab=xlab,ylab=ylab,pch=16)
    lines(x[data[,1]=="A"],data[,2][data[,1]=="A"])
    lines(x[data[,1]=="B"],data[,2][data[,1]=="B"])
    lines(c(sum(data[,1]=="A")+0.5,sum(data[,1]=="A")+0.5),c(min(data[,2])-5,max(data[,2])+5),lty=2)
    mtext("A",side=3,at=(sum(data[,1]=="A")+1)/2)
    mtext("B",side=3,at=(sum(data[,1]=="A")+(sum(data[,1]=="B")+1)/2))
  }
  
  
  if(design=="ABA"){
    plot(x,data[,2],xlab=xlab,ylab=ylab,pch=16)
    lines(x[data[,1]=="A1"],data[,2][data[,1]=="A1"])
    lines(x[data[,1]=="B1"],data[,2][data[,1]=="B1"])
    lines(x[data[,1]=="A2"],data[,2][data[,1]=="A2"])
    lines(c(sum(data[,1]=="A1")+0.5,sum(data[,1]=="A1")+0.5),c(min(data[,2])-5,max(data[,2])+5),lty=2)
    lines(c(sum(data[,1]=="A1")+sum(data[,1]=="B1")+0.5,sum(data[,1]=="A1")+sum(data[,1]=="B1")+0.5),c(min(data[,2])-5,max(data[,2])+5),lty=2)
    mtext("A",side=3,at=(sum(data[,1]=="A1")+1)/2)
    mtext("B",side=3,at=(sum(data[,1]=="A1")+(sum(data[,1]=="B1")+1)/2))
    mtext("A",side=3,at=(sum(data[,1]=="A1")+sum(data[,1]=="B1")+(sum(data[,1]=="A2")+1)/2))
  }
  
  if(design=="ABAB"){
    plot(x,data[,2],xlab=xlab,ylab=ylab,pch=16)
    lines(x[data[,1]=="A1"],data[,2][data[,1]=="A1"])
    lines(x[data[,1]=="B1"],data[,2][data[,1]=="B1"])
    lines(x[data[,1]=="A2"],data[,2][data[,1]=="A2"])
    lines(x[data[,1]=="B2"],data[,2][data[,1]=="B2"])
    lines(c(sum(data[,1]=="A1")+0.5,sum(data[,1]=="A1")+0.5),c(min(data[,2])-5,max(data[,2])+5),lty=2)
    lines(c(sum(data[,1]=="A1")+sum(data[,1]=="B1")+0.5,sum(data[,1]=="A1")+sum(data[,1]=="B1")+0.5),c(min(data[,2])-5,max(data[,2])+5),lty=2)
    lines(c(sum(data[,1]=="A1")+sum(data[,1]=="B1")+sum(data[,1]=="A2")+0.5,sum(data[,1]=="A1")+sum(data[,1]=="B1")+sum(data[,1]=="A2")+0.5),c(min(data[,2])-5,max(data[,2])+5),lty=2)
    mtext("A",side=3,at=(sum(data[,1]=="A1")+1)/2)
    mtext("B",side=3,at=(sum(data[,1]=="A1")+(sum(data[,1]=="B1")+1)/2))
    mtext("A",side=3,at=(sum(data[,1]=="A1")+sum(data[,1]=="B1")+(sum(data[,1]=="A2")+1)/2))
    mtext("B",side=3,at=(sum(data[,1]=="A1")+sum(data[,1]=="B1")+sum(data[,1]=="A2")+(sum(data[,1]=="B2")+1)/2))
  }
  
  if(design=="MBD"){
    N<-ncol(data)/2
    par(mfrow=c(N,1))
    for(it in 1:N){
      plot(x,data[,it*2],xlab="",ylab=ylab,pch=16)
      lines(x[data[,(it*2)-1]=="A"],data[,it*2][data[,(it*2)-1]=="A"])
      lines(x[data[,(it*2)-1]=="B"],data[,it*2][data[,(it*2)-1]=="B"])
      lines(c(sum(data[,(it*2)-1]=="A")+0.5,sum(data[,(it*2)-1]=="A")+0.5),c(min(data[,it*2])-5,max(data[,it*2])+5),lty=2)
      mtext("A",side=3,at=(sum(data[,(it*2)-1]=="A")+1)/2)
      mtext("B",side=3,at=(sum(data[,(it*2)-1]=="A")+(sum(data[,(it*2)-1]=="B")+1)/2))
    }
    title(xlab=xlab,pch=16)
    
    par(mfrow=c(1,1))
  }


}
