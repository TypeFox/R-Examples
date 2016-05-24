plot.TPidm <-
function(x,chosen.tr="ALL", ...){
    if(!inherits(x,"TPidm")) stop("'x' must be of class 'TPidm'.")
    times<-x$times
    s<-x$s
    t<-x$t
    CI<-x$CI
    
    if(chosen.tr[1]=="ALL") chosen.tr<- x$p.trans
    
    
    
    # CI==TRUE
    if(CI==TRUE){
      all.probs<-x$all.probs[,,chosen.tr]
      
      # option 1: individual curves with CI
      if(length(chosen.tr)==1){
        tit<-paste("p",chosen.tr)
        plot(times,all.probs[,"probs"],type="s",xlab="Time",ylab="Probability",ylim=c(0,1),main=tit, ...)
        lines(times,all.probs[,"lower"],type="l",lty=3,col="red")
        lines(times,all.probs[,"upper"],type="l",lty=3,col="red")
        
      }else if(length(chosen.tr)>1 & length(chosen.tr)<=3){
        N<-seq(1,length(chosen.tr),1)
        tpg<-layout(matrix(N,1,length(chosen.tr),byrow=TRUE))
        layout.show(tpg)
        
        for(i in 1:length(chosen.tr)){
          tit<-paste("p",chosen.tr[[i]])
          plot(times,all.probs[,"probs",chosen.tr[[i]]],type="s",xlab="Time",ylab="Probability",ylim=c(0,1),main=tit, ...)
          lines(times,all.probs[,"lower",chosen.tr[[i]]],type="l",lty=3,col="red")
          lines(times,all.probs[,"upper",chosen.tr[[i]]],type="l",lty=3,col="red") 
        }
      } else {
        
        n.cols<-3
        n.rows<-ceiling(length(chosen.tr)/n.cols)
        N<-seq(1,length(chosen.tr),1)
        if(length(N)<n.cols*n.rows){
          m<-n.cols*n.rows-length(N)
          N<-c(N,rep(0,m))
        }
        tpg<-layout(matrix(N,n.rows,n.cols,byrow=TRUE))
        layout.show(tpg)
        
        for(i in 1:length(chosen.tr)){
          tit<-paste("p",chosen.tr[[i]])
          plot(times,all.probs[,"probs",chosen.tr[[i]]],type="s",xlab="Time",ylab="Probability",ylim=c(0,1),main=tit, ...)
          lines(times,all.probs[,"lower",chosen.tr[[i]]],type="l",lty=3,col="red")
          lines(times,all.probs[,"upper",chosen.tr[[i]]],type="l",lty=3,col="red") 
        }
      }
      
      if(s==0){
        title("State Occupation Probabilities", outer=TRUE,  line=-1)
      } else{
        title("Transition Probabilities", outer=TRUE,  line=-1)
      }
      
    }else{
      
      
      all.probs<-x$all.probs[,,chosen.tr]
      
      # option 2: individual curves without CI
      
      if(length(chosen.tr)==1){
        tit<-paste("p",chosen.tr)
        plot(times,all.probs,type="s",xlab="Time",ylab="Probability",ylim=c(0,1),main=tit, ...)
        
      }else if(length(chosen.tr)>1 & length(chosen.tr)<=3){
        N<-seq(1,length(chosen.tr),1)
        tpg<-layout(matrix(N,1,length(chosen.tr),byrow=TRUE))
        layout.show(tpg)
        
        for(i in 1:length(chosen.tr)){
          tit<-paste("p",chosen.tr[[i]])
          plot(times,all.probs[,chosen.tr[[i]]],type="s",xlab="Time",ylab="Probability",ylim=c(0,1),main=tit, ...)
          
        }
      } else {
        
        n.cols<-3
        n.rows<-ceiling(length(chosen.tr)/n.cols)
        N<-seq(1,length(chosen.tr),1)
        if(length(N)<n.cols*n.rows){
          m<-n.cols*n.rows-length(N)
          N<-c(N,rep(0,m))
        }
        tpg<-layout(matrix(N,n.rows,n.cols,byrow=TRUE))
        layout.show(tpg)
        
        for(i in 1:length(chosen.tr)){
          tit<-paste("p",chosen.tr[[i]])
          plot(times,all.probs[,chosen.tr[[i]]],type="s",xlab="Time",ylab="Probability",ylim=c(0,1),main=tit, ...)
          
        }
      }
      if(s==0){
        title("State Occupation Probabilities", outer=TRUE,  line=-1)
      } else{
        title("Transition Probabilities", outer=TRUE,  line=-1)
      }
      
    } # if CI==TRUE
  }
