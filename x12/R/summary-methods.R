setMethod("summary",
    signature(object = "x12Output"),
    function (object, fullSummary=FALSE, spectra.detail=FALSE, almostout=FALSE, rsd.autocorr=NULL, quality.stat=FALSE, likelihood.stat=FALSE, aape=FALSE, id.rsdseas=FALSE, slidingspans=FALSE, history=FALSE, identify=FALSE, print=TRUE) 
    {
      if(length(object@dg)>0){
        if(print)
          summaryworkhorse(object@dg,fullSummary=fullSummary,spectra.detail=spectra.detail,
              almostout=almostout,rsd.autocorr=rsd.autocorr,quality.stat=quality.stat,
              likelihood.stat=likelihood.stat,aape=aape,id.rsdseas=id.rsdseas,
              slidingspans=slidingspans,history=history,identify=identify)
        invisible(summary.output.workhorse(object@dg,fullSummary=fullSummary,
                spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,
                quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,
                id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,identify=identify))
      }else{
        cat("You need to run x12 before viewing a summary!\n")
      }
    }
)
setMethod("summary",
    signature(object = "x12Single"),
    function (object, fullSummary=FALSE, spectra.detail=FALSE, almostout=FALSE, rsd.autocorr=NULL,
        quality.stat=FALSE, likelihood.stat=FALSE, aape=FALSE, id.rsdseas=FALSE,slidingspans=FALSE,
        history=FALSE, identify=FALSE, oldOutput=NULL,print=TRUE) 
    { 
      ############
#		sumout<-unlist(summary.output[,1])
#		names.sumout<-unique(sumout)
#		length.names.sumout<-unlist(lapply(names.sumout,function(x)length(grep(x,sumout))))
      ###########
#	summary.output<-summary(object@x12Output,fullSummary=fullSummary,spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,identify=identify)
      if(is.null(oldOutput)){
        if(!is.null(object@tsName)){  
          if(print){
            cat("--------------------------  ",object@tsName,"  ------------------------------------\n")
            cat("-----------------------------------------------------------------------------------\n")
          }
          summary.output<-summary(object@x12Output,fullSummary=fullSummary,
              spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,
              quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,
              id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,
              identify=identify,print=print)
          names(summary.output)[2]<-object@tsName
        }else{  
          if(print){
            cat("--------------------------  Rout  ------------------------------------\n")
            cat("-----------------------------------------------------------------------------------\n")
          }
          summary.output<-summary(object@x12Output,fullSummary=fullSummary,
              spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,
              quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,
              id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,
              identify=identify,print=print)
        }
      }else{
        nprev <- min(length(object@x12OldOutput),oldOutput)
        if(!is.null(object@tsName)){
          if(print){
            cat("--------------------------  ",object@tsName,"  ------------------------------------\n")
            cat("-----------------------------------------------------------------------------------\n")
          }
          summary.output<-summary(object@x12Output,fullSummary=fullSummary,
              spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,
              quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,
              id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,
              identify=identify,print=print)
          names(summary.output)[2]<-object@tsName
        }else{
          if(print){
            cat("--------------------------  Rout  ------------------------------------\n")
            cat("-----------------------------------------------------------------------------------\n")
          }
          summary.output<-summary(object@x12Output,fullSummary=fullSummary,
              spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,
              quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,
              id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,
              identify=identify,print=print)
        }
        if(nprev>0){
          for(i in nprev:1){
            if(i==length(object@x12OldOutput))
              TF <- !identical(object@x12Output,object@x12OldOutput[[i]])
            else if(i!=nprev)
              TF <- !identical(object@x12OldOutput[[i]],object@x12OldOutput[[i+1]])
            else
              TF <- TRUE
            if(TF){
              if(print)
                cat("\n---------------------------  RUN  ",i,"  ----------------------------------------\n")
              summary.oldout<-summary(object@x12OldOutput[[i]],fullSummary=fullSummary,
                  spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,
                  quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,
                  id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,
                  identify=identify,print=print)
              names(summary.oldout)<-names(summary.output)
              summary.oldout<-rbind(c(paste("OLD OUTPUT",i),paste("Run",i)),summary.oldout)
              summary.output<-rbind(summary.output,summary.oldout)
            }else{
              if(print)
                cat("--- No valid previous runs. ---\n")
            }
          }  
        }
      }
      invisible(summary.output)
    }
)
setMethod("summary",
    signature(object = "x12Batch"),
    function (object, fullSummary=FALSE, spectra.detail=FALSE, almostout=FALSE, rsd.autocorr=NULL,
        quality.stat=FALSE, likelihood.stat=FALSE, aape=FALSE, id.rsdseas=FALSE,
        slidingspans=FALSE,history=FALSE,identify=FALSE,oldOutput=NULL,print=TRUE) 
    {
      summary.output<-list()
      for(i in 1:length(object@x12List)){
        if(print)
          cat("-----------------------------------------------------------------------------------\n")
        summary.output[[length(summary.output)+1]]<-summary(object@x12List[[i]],fullSummary=fullSummary,spectra.detail=spectra.detail,almostout=almostout,rsd.autocorr=rsd.autocorr,quality.stat=quality.stat,likelihood.stat=likelihood.stat,aape=aape,id.rsdseas=id.rsdseas,slidingspans=slidingspans,history=history,identify=identify,oldOutput=oldOutput,print=print)
      }
      if(is.null(oldOutput)){
        summary.output<-fun.table(summary.output,object=object)		
        invisible(summary.output)
      }else{
        new.out<-list()
        old.out<-list()
        for(i in 1:length(object@x12List)){
#		cat("i=",i,"\n")
          if(any(grep("OLD OUTPUT",summary.output[[i]][,1]))){		
            new.out[[i]]<-summary.output[[i]][1:(min(grep("OLD OUTPUT",summary.output[[i]][,1]))-1),]
            old.out.sub<-list()
            ind.old.out<-grep("OLD OUTPUT",summary.output[[i]][,1])
            if(length(ind.old.out)>1){
              for(j in 1:(length(ind.old.out)-1)){
                old.out.sub[[length(old.out.sub)+1]]<- summary.output[[i]][ind.old.out[j]:(ind.old.out[j+1]-1),]		
              }}
            old.out.sub[[length(old.out.sub)+1]]<- summary.output[[i]][ind.old.out[length(ind.old.out)]:length(summary.output[[i]][,1]),]		
            old.out[[length(old.out)+1]]<-old.out.sub
          }else{
            new.out[[i]]<-summary.output[[i]]	
            old.out[[length(old.out)+1]]<-NA
          }
        }  #
        summary.output<-fun.table(new.out,object=object)
        noo.values<-sort(unique(grep("OLD OUTPUT",unlist(old.out),value=TRUE)),decreasing=TRUE)
        noo <- length(noo.values)
        new.list<-list()
        new.list<-lapply(1:length(old.out),function(x)
              new.list[[length(new.list)+1]]<-rep(data.frame(NA),noo))
        
        for(i in 1:length(old.out)){
          if(noo>0){
            for(j in 1:noo){
              if(any(grepl(noo.values[j],old.out[[i]])))	
                new.list[[i]][[j]]<-old.out[[i]][[grep(noo.values[j],old.out[[i]])]]
            }
          }
        }	
        old.out<-new.list
        
        new.list<-list()
        if(noo>0){
          for(i in 1:noo){
            sub.new.list<-list()	
            new.list[[length(new.list)+1]]<-lapply(old.out,function(x){
                  if(all(is.na(x[[i]]))){
                    x[[i]]<-list(c(paste("OLD OUTPUT",c(noo:1)[i]),"No previous run"))
                    x[[i]]<-unlist(x[[i]])
                  }
                  sub.new.list[[length(sub.new.list)+1]]<-x[[i]]	
                })	
          }
          
          for(i in 1:noo){
            summary.output<-rbind(summary.output,fun.table(new.list[[i]],object=object))	
          }
        }
      }
      invisible(summary.output)
    }
)


fun.table <- function(y,object){		
  row.names<-unique(unlist(lapply(y,function(x){
                if(is.character(x))
                  x[1]
                else
                  x[,1]})))
  new.out<-cbind(row.names,rep(NA,length(row.names)))
  col.names<-vector()
  for(i in 1:length(object@x12List)){
    new.out2<-new.out
    if(is.character(y[[i]])){
      new.out2<-new.out2[,2]
      col.names[length(col.names)+1]<-object@x12List[[i]]@tsName	
    }
    else{
      new.out2<-apply(new.out2,1,function(x){
            if(x[1]%in%y[[i]][,1])								
              x[2]<-y[[i]][which(y[[i]][,1]%in%x[1]),2]	
            else
              x[2]<-x[2]
          })
      col.names[length(col.names)+1]<-names(y[[i]])[2]
    }
    if(i==1)
      new.out3<-data.frame(cbind(new.out[,1],unlist(new.out2)))
    else
      new.out3<-cbind(new.out3,unlist(new.out2))
  }
  y<-as.data.frame(new.out3)
  colnames(y)<-c("DIAGNOSTICS",col.names)
  g <- grep("variable",as.character(y[,1]))
  if(length(g)>1){
	  ind1 <- 1:(g[1]-1)
	  y2 <- y[ind1,]
	  ind2 <- g
	  y <- rbind(y2,y[ind2,],y[-c(ind1,ind2),])
	  rownames(y) <- 1:nrow(y)
  }
  
  return(y)
}
