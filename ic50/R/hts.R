hts.96<-function(indir=".",plates=2,
                  measure=NULL,control=NULL,dilution=NULL,
                  inhib=NULL,normalize="mean",graphics="mean",outdir="./results"){
  files<-list.files(path=indir)
  if(!file.exists(paste(outdir,"/",sep=""))){
    warning(paste("New output directory",outdir,"created."))
    dir.create(outdir)
  }
  pdf(file=paste(outdir,"/dose_response_curves.pdf",sep=""),height=8,width=11)
  for(i in seq(1,length(files)-plates+1,plates)){
    cat("Evaluating data in",files[i:(i+plates-1)],"\n")
    results_cp<-ic50.96(files=paste(indir,"/",files[i:(i+plates-1)],sep=""),
                        measure=measure,control=control,dilution=dilution,
                        inhib=inhib,normalize=normalize,graphics=graphics,outdir=NULL)
    filenames<-rep(files[i],length(results_cp$ic50))
    if(i==1) results<-cbind(filenames,results_cp)
    else results<-rbind(results,cbind(filenames,results_cp))
  }
  dev.off()
  names(results)<-c("first_file",names(results_cp))
  if(!is.null(outdir)) write.table(results,file=paste(outdir,"/ic50.txt",sep=""),sep="\t",row.names=FALSE)
  cat("Results written to folder ",outdir,". Thank you.\n\n",sep="")
  return(results)
}

hts.384<-function(indir=".",plates=2,
                  measure=NULL,control=NULL,dilution=NULL,
                  inhib=NULL,normalize="single",graphics="mean",outdir="./results"){
  files<-list.files(path=indir)
  if(!file.exists(paste(outdir,"/",sep=""))){
    warning(paste("New output directory",outdir,"created."))
    dir.create(outdir)
  }
  pdf(file=paste(outdir,"/dose_response_curves.pdf",sep=""),height=8,width=11)
  for(i in seq(1,length(files)-plates+1,plates)){
    cat("Evaluating data in",files[i:(i+plates-1)],"\n")
    results_cp<-ic50.384(files=paste(indir,"/",files[i:(i+plates-1)],sep=""),
                         measure=measure,control=control,dilution=dilution,
                         inhib=inhib,normalize=normalize,graphics=graphics,outdir=NULL)
    filenames<-rep(files[i],length(results_cp$ic50))
    if(i==1) results<-cbind(filenames,results_cp)
    else results<-rbind(results,cbind(filenames,results_cp))
  }
  dev.off()
  names(results)<-c("first_file",names(results_cp))
  if(!is.null(outdir)) write.table(results,file=paste(outdir,"/ic50.txt",sep=""),sep="\t",row.names=FALSE)
  cat("Results written to folder ",outdir,". Thank you.\n\n",sep="")
  return(results)
}
