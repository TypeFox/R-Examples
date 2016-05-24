#The workhorse
evaluation<-function(y,dilution,inhib,outdir,file,graphics){
  cpnames<-names(y)
  ncp<-length(cpnames)
  ic50mean<-ygr_mean<-cl<-cu<-area<-maxsd<-cv<-ncon<-nrep<-numeric(0)
  for(i in 1:ncp){
    nrep[i]<-length(y[[i]][,1])
    ncon[i]<-length(dilution[[i]])
  }
  stddev<-ic50<-list()

  
  #Create a grid
  x<-list()
  y_mean<-numeric(0)
  for(cp in 1:ncp) x[[cp]]<-log10(dilution[[cp]])
  xgr<-ygr<-list()
  for(cp in 1:ncp){
    ord<-order(x[[cp]])
    x[[cp]]<-x[[cp]][ord]
    xgr[[cp]]<-linearCurve(x[[cp]],y[[cp]][1,])[[1]]
    
    ygr[[cp]]<-matrix(nrow=nrep[cp],ncol=length(xgr[[cp]]))
    for(rep in 1:nrep[cp]){
      y[[cp]][rep,]<-y[[cp]][rep,ord]
      ygr[[cp]][rep,]<-linearCurve(x[[cp]],y[[cp]][rep,])[[2]]
    }
  }
  
  #Evaluation
  for(cp in 1:ncp){
    stddev[[cp]]<-ic50[[cp]]<-numeric(0)

    for(rep in 1:nrep[cp]){
      if(min(ygr[[cp]][rep,])<=1-inhib[cp] && ygr[[cp]][rep,1]>1-(.5*inhib[cp])) ic50[[cp]][rep]<-preimage(1-inhib[cp],xgr[[cp]],ygr[[cp]][rep,])
      else ic50[[cp]][rep]<-NA #IC50s from single curves
    }
    nval<-sum(!is.na(ic50[[cp]]))
    
    for(cn in 1:length(ygr[[cp]][1,])) ygr_mean[cn]<-mean(ygr[[cp]][,cn])    
    if(min(ygr_mean)<=1-inhib[cp] && ygr_mean[1]>1-(.5*inhib[cp])) ic50mean[cp]<-preimage(1-inhib[cp],xgr[[cp]],ygr_mean)
    else ic50mean[cp]<-NA #IC50s
 
    for(i in 1:ncon[cp]) stddev[[cp]][i]<-sd(y[[cp]][,i])
    maxsd[cp]<-max(stddev[[cp]]) #Maximum standard deviation
    if(2*nval>=length(ic50[[cp]])){
      cv[cp]<-sd(ic50[[cp]],na.rm=TRUE)/abs(ic50mean[cp]) #Coeff of variation
    }
    else cv[cp]<-NA

    if(2*nval>=length(ic50[[cp]]) && var(ic50[[cp]],na.rm=TRUE)!=0 && sum(!is.na(ic50[[cp]]))>=2){
      cl[cp]<-ic50mean[cp]-qt(p=.975,df=nval)*sqrt(var(ic50[[cp]],na.rm=TRUE)/nval) #Confidence interval
      cu[cp]<-ic50mean[cp]+qt(p=.975,df=nval)*sqrt(var(ic50[[cp]],na.rm=TRUE)/nval)
    }
    else cl[cp]<-cu[cp]<-NA
 
  }
  
  results<-data.frame(cpnames,round(10^ic50mean,4),round(10^cl,4),round(10^cu,4),round(maxsd,4),round(cv,4))
  colnames(results)<-c("compound","ic50","clow","cup","maxsd","cv")

  
  #Create graphical output
  if(!is.null(outdir)){
    if(!file.exists(paste(outdir,"/",sep=""))){
      warning(paste("New output directory",outdir,"created."))
      dir.create(outdir)
    }
    pdf(file=paste(outdir,"/dose_response_curves.pdf",sep=""),height=8,width=11)
  }
  for(cp in 1:ncp){
    if(graphics=="mean") meanplot(x[[cp]],y[[cp]],ic50mean[cp],cl[cp],cu[cp],stddev[[cp]],file,cpnames[cp])
    if(graphics=="fitted") fittedplot(x[[cp]],y[[cp]],ic50mean[cp],cl[cp],cu[cp],file,cpnames[cp])
    if(graphics=="single") singleplot(x[[cp]],y[[cp]],ic50mean[cp],cl[cp],cu[cp],file,cpnames[cp])
  }

  
  #Write results and data to files
  if(!is.null(outdir)){
    dev.off()
    write.table(results,file=paste(outdir,"/ic50.txt",sep=""),sep="\t",row.names=FALSE)
    file.remove(paste(outdir,"/measurement.txt",sep=""))
    for(cp in 1:ncp){
      cat(paste(cpnames[cp],"\n"),file=paste(outdir,"/measurement.txt",sep=""),append=TRUE)
      write.table(round(y[[cp]],4),file=paste(outdir,"/measurement.txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
      cat("\n",file=paste(outdir,"/measurement.txt",sep=""),append=TRUE)
    }
    cat("Results written to folder ",outdir,". Thank you.\n\n",sep="")
  }
  return(results)
}
