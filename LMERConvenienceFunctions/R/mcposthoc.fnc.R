mcposthoc.fnc<-function(model,var,two.tailed=TRUE,
    mcmc=FALSE,nsim=10000,ndigits=4,mc.cores=1,
    verbosity=1,...){
  # check if parallel is available
  if(!"parallel" %in% .packages(all.available = TRUE)){
    stop("package parallel not available on this machine\n")
  }

  if(mcmc){
    stop("This argument is depricated. Please see \"http://stackoverflow.com/questions/19199713/lme4-and-languager-compatibility-error-input-model-is-not-a-mer-object\" for possible avenues to get p-values.\n")
  }

  if(!is.list(var)){
    stop("argument \"var\" must be a list\n")
  }
  
  if(is.null(names(var))){
    stop("\"var\" must be a named list\n")
  }

  data<-model@frame
  
  # create list of posthoc tests to do
  for(k in 1:length(names(var))){
    x<-t(unique(eval(parse(text=paste("data[,var[[k]]]",collapse=",")))))
    colnames(x)<-paste(names(var)[k],1:ncol(x),sep="__")
                if(length(var[[k]])==1)rownames(x)<-var[k]
    for(kk in 1:length(x[1,])){
      if(kk==1){
        if(k==1){
          ph.list<-list()
        }
        tmp<-x[,kk]
        names(tmp)<-rownames(x)
        ph.list[[colnames(x)[kk]]]<-tmp
      }else{
        tmp<-x[,kk]
        names(tmp)<-rownames(x)
        ph.list[[colnames(x)[kk]]]<-tmp
      }

    }
  }

  # create vector of collapsed names for cat
  my.nms<-vector("character")
  for(ff in 1:length(ph.list)){
    my.nms<-c(my.nms,paste("\"",gsub(" ","",paste(names(ph.list[[ff]]),ph.list[[ff]],collapse="_")),"\"",sep=""))
  }

  if(!mcmc){
    my.update<-function(model,ph.list.element,verbosity){
    # do posthoc for an element of ph.list --> relevel and update model
    for(j in 1:(length(names(ph.list.element)))){
      predictor<-names(ph.list.element)[j]
      data[,predictor]<-relevel(data[,predictor], as.character(ph.list.element[j]))
    }
    myjob<-paste("\"",gsub(" ","",paste(names(ph.list.element),
    ph.list.element,collapse="_")),"\"",sep="")
    oneofwhat<-paste("(",which(my.nms==myjob)," of ",length(my.nms),")",sep="")
    if(verbosity>0)cat("processing job",myjob,oneofwhat,"...\n")
      return(invisible(summary(update(model,.~.,data=data))$coefficients))
    }
  
    ph.smrys<-mclapply(ph.list,FUN=function(x)my.update(model=model,
    ph.list.element=x,verbosity=verbosity),mc.cores=mc.cores,...)

    # create ph.names
    ph.names<-gsub("(.*)__.*","\\1",names(ph.smrys))
    for(ii in 1:length(ph.names)){
      ph.names[ii]<-paste(ph.names[ii],paste(ph.list[[ii]],collapse="_"),sep="_")
    }
    names(ph.smrys)<-ph.names
  
    for(mm in 1:length(ph.names)){
      # store summary data in object temp
      temp<-as.data.frame(ph.smrys[[mm]])
  
      # get rank of X matrix
      rank.X=qr(getME(model,"X"))$rank
  
      # get upper-bound df
      udf<-nrow(model@frame)-rank.X
  
      # get lower-bound df
      model.ranef<-ranef(model)
      lower.bound<-0
      for(nn in 1:length(names(model.ranef))){
        dims<-dim(model.ranef[[nn]])
        lower.bound<-lower.bound+dims[1]*dims[2]
      }
      ldf<-nrow(model@frame)-rank.X-lower.bound
         
      multip<-ifelse(two.tailed,2,1)
      temp[,"udf"]<-udf
      temp[,"ldf"]<-ldf
      if(as.vector(model@call[1])=="glmer()"){
        temp[,"upper.p.val"]<-round(multip*(1-pt(abs(temp[,"z value"]),udf)),ndigits)
        temp[,"lower.p.val"]<-round(multip*(1-pt(abs(temp[,"z value"]),ldf)),ndigits)
      }else{
        temp[,"upper.p.val"]<-round(multip*(1-pt(abs(temp[,"t value"]),udf)),ndigits)
        temp[,"lower.p.val"]<-round(multip*(1-pt(abs(temp[,"t value"]),ldf)),ndigits)
      }
        
      # put temp back in summary object
      ph.smrys[[mm]]<-temp
    }
  }

  if(mcmc){
	### DEPRICATED SINCE VERSION 2.3
    my.mc.update<-function(model,ph.list.element,nsim=nsim,ndigits=ndigits,verbosity=verbosity,...){
      # do posthoc for an element of ph.list --> relevel and update model
      for(j in 1:(length(names(ph.list.element)))){
        predictor<-names(ph.list.element)[j]
        data[,predictor]<-relevel(data[,predictor],as.character(ph.list.element[j]))
      }
      myjob<-paste("\"",gsub(" ","",paste(names(ph.list.element), ph.list.element,collapse="_")),"\"",sep="")
      oneofwhat<-paste("(",which(my.nms==myjob)," of ",length(my.nms),")",sep="")
      if(verbosity>0)cat("processing job",myjob,oneofwhat,"...\n")
      updated.model<-update(model,.~.,data=data)
	  pvals.fnc<-NULL
      model.mcmc<-pvals.fnc(object=updated.model,nsim=nsim, ndigits=ndigits,withMCMC=FALSE,addPlot=FALSE,...)
      return(invisible(model.mcmc$fixed))
    }
  
    ph.smrys<-mclapply(ph.list,FUN=function(x)my.mc.update(model=model,ph.list.element=x,nsim=nsim,ndigits=ndigits,verbosity=verbosity),mc.cores=mc.cores,...)
  
    # create ph.names
    ph.names<-gsub("(.*)__.*","\\1",names(ph.smrys))
    for(ii in 1:length(ph.names)){
      ph.names[ii]<-paste(ph.names[ii],paste(ph.list[[ii]],collapse="_"),sep="_")
    }
    names(ph.smrys)<-ph.names
  }

  res<-list(n=nrow(model@frame),var=var,summaries=ph.smrys)
  class(res)<-"mcposthoc"
  # return results
  return(res)
}

