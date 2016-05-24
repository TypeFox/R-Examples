

getPSS<-function(riskmeta,plot=TRUE,binsize=400){

      
      nodeseq<-Getenvir("getnodeseq")
      strnew<-Getenvir("getStr")
      
      riskmetainnet<-intersect(nodeseq,riskmeta)
      riskmetainnet<-unique(riskmetainnet)
      
      
      
      infor<-list()
      for(i in 1:length(strnew[,1]))
         {

            infor_i<-list(nodeId="000",meanvalue=0,diffnode=0)
     
            infor_i$nodeId<-strnew[i,1]
            infor_i$meanvalue<-strnew[i,2]
    
            if(strnew[i,1]%in%riskmetainnet)
              {
                 infor_i$diffnode<-1
              }
     
     
            infor[[i]]<-infor_i
   
          }
      nodeId<-sapply(infor,function(x) x$nodeId)
      meanvalue<-sapply(infor,function(x) x$meanvalue)
      diffnode<-sapply(infor,function(x) x$diffnode)
      
      meanstr<-data.frame(nodeId=nodeId,meanvalue=meanvalue,diffnode=diffnode)
      meanstr<-meanstr[order(meanstr[,2],decreasing = FALSE),]
      
      pss<-performpcls(meanstr[,"meanvalue"],meanstr[,"diffnode"])
      CGNB<-1-pss
      re_out<-data.frame(riskmeta=meanstr[,"diffnode"],meanstrvalue=meanstr[,"meanvalue"],pss=pss,CGNB=CGNB,stringsAsFactors=FALSE)
      rownames(re_out)=meanstr[,"nodeId"]
     

     
     
     
     
      if(plot)
       {
         splitter=ceiling(1:length(re_out$riskmeta)/binsize)
			   de=sapply(split(re_out$riskmeta,splitter),mean)
			   binlen=sapply(split(as.numeric(re_out$meanstrvalue),splitter),mean)
		     plot(binlen,de,xlab=paste("Bias Data in ",binsize," metabolite bins",sep=""),ylab="Proportion risk metabolite")
		     lines(re_out$meanstrvalue,re_out$pss,col=3,lwd=2)
		
		   }#####plot
 
     
      return(re_out)

      }###########function()
