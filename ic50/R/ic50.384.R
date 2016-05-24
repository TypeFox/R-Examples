ic50.384<-function(files,measure=NULL,control=NULL,dilution=NULL,inhib=NULL,normalize="single",graphics="mean",outdir="./results"){
  defaultfiles.write()
  data<-data_merged<-dillist<-list()
  nconc<-ncontr<-numeric(0)

  exc<-catch.exceptions.384(files,measure,control,dilution,inhib,graphics,normalize)
  measure<-exc[[1]]
  control<-exc[[2]]
  dilution<-exc[[3]]
  
  measnames<-measure[,1]
  cpnames<-levels(as.factor(measnames))
  
  for(i in 1:length(measnames)) nconc[i]<-sum(as.character(dilution[i,-1])!="")
  for(i in 1:length(measnames)) ncontr[i]<-sum(as.character(control[i,-1])!="")
  for(j in 1:length(measnames)){for(i in 1:length(cpnames)){
    if(measnames[j]==cpnames[i]){
      dillist[[i]]<-as.numeric(as.character(dilution[j,-1]))
    }
  }}

  for(file in 1:length(files)){
    rawdata<-read.delim(files[file],sep="\t",header=FALSE)
    if(any(dim(rawdata)!=c(16,24))) stop("384 (16 x 24) wells needed.")
    for(ms in 1:length(measnames)){
      controlrow<-thisrow<-numeric(0)
      if(normalize=="single"){
        if(length(measure[ms,])!=length(control[ms,])) stop("Need one control well per measure well.")
        for(cn in 1:nconc[ms]) thisrow<-c(thisrow,eval(parse(text=paste("rawdata[",measure[ms,cn+1],"]/rawdata[",control[ms,cn+1],"]",sep=""))))
      }
      if(normalize=="mean"){
        for(cn in 1:ncontr[ms]) controlrow<-c(controlrow,eval(parse(text=paste("rawdata[",control[ms,cn+1],"]",sep=""))))
        if(mean(controlrow)==1){
          thisrow<-rep(0,nconc[ms])
        }
        else{
          for(cn in 1:nconc[ms]) thisrow<-c(thisrow,eval(parse(text=paste("rawdata[",measure[ms,cn+1],"]/mean(controlrow)",sep=""))))
        }
      }
      data[[(file-1)*length(measnames)+ms]]<-thisrow
      names(data)[(file-1)*length(measnames)+ms]<-measnames[ms]
    }
  }
  
  for(cp in 1:length(cpnames)) data_merged[[cp]]<-numeric(0)
  for(ms in 1:length(data)){for(cp in 1:length(cpnames)){
    if(names(data)[ms]==cpnames[cp]){
      data_merged[[cp]]<-rbind(data_merged[[cp]],data[[ms]])
    }
  }}
  if(is.null(inhib)) inhib<-rep(0.5,length(cpnames))
  names(dillist)<-names(data_merged)<-cpnames

  evaluation(data_merged,dillist,inhib,outdir,files[1],graphics)
}



catch.exceptions.384<-function(files,measure,control,dilution,inhib,graphics,normalize){
  if(is.null(measure)){
    warning("No measure wells specified, .last384_measure.txt is used",call.=FALSE)
    measure<-as.matrix(read.delim(".last384_measure.txt",colClasses="character",header=FALSE))
  }
  else{
    measure<-as.matrix(read.delim(measure,colClasses="character",header=FALSE))
    write.table(measure,file=".last384_measure.txt",row.names=FALSE,col.names=FALSE,sep="\t")
  }
  if(length(levels(as.factor(as.character(measure[,-1]))))!=length(as.character(measure[,-1]))){
    stop("Multiple use of wells.")
  }
  if(is.null(control)){
    warning("No control wells specified, .last384_control.txt is used",call.=FALSE)
    control<-as.matrix(read.delim(".last384_control.txt",colClasses="character",header=FALSE))
  }
  else{
    control<-as.matrix(read.delim(control,colClasses="character",header=FALSE))
    write.table(control,file=".last384_control.txt",row.names=FALSE,col.names=FALSE,sep="\t")
  }
  if(is.null(dilution)){
    warning("No dilutions specified, .last384_dilution.txt is used",call.=FALSE)
    dilution<-as.matrix(read.delim(".last384_dilution.txt",colClasses="character",header=FALSE))
  }
  else{
    dilution<-as.matrix(read.delim(dilution,colClasses="character",header=FALSE))
    write.table(dilution,file=".last384_dilution.txt",row.names=FALSE,col.names=FALSE,sep="\t")
  }
  if(any(is.na(as.numeric(dilution[,-1])))) stop("Numeric concentrations needed in dilution file.")

  measnames<-measure[,1]
  cpnames<-levels(as.factor(measnames))
  if(any(measnames!=control[,1])) stop("Row names differ for measure wells and controls.")
  if(any(measnames!=dilution[,1])) stop("Row names differ for measure wells and dilutions.")

  nconc<-numeric(0)
  for(i in 1:length(measnames)){
    nconc[i]<-sum(as.character(dilution[i,-1])!="")
    if(sum(as.character(measure[i,-1])!="")!=nconc[i]){
      stop(paste("Need one well per concentration for compound ",measnames[i],".",sep=""))
    }
  }
  for(i in 1:(length(measnames)-1)){for(j in (i+1):length(measnames)){
    if(measnames[i]==measnames[j]){
      if(nconc[i]!=nconc[j]){
        stop(paste("Dilution rows differ for compound ",measnames[i],".",sep=""))
      }
      if(any(as.numeric(as.character(dilution[i,-1])) != as.numeric(as.character(dilution[j,-1])),na.rm=TRUE)){
        cat(dilution[i,-1]," ",dilution[j,-1])
        stop(paste("Dilution rows differ for compound ",measnames[i],".",sep=""))
      }
    }
  }}
  if(!is.null(inhib)){
    if(length(inhib)!=length(cpnames)) stop("Need one inhibitory percentage per compound.")
    if(max(inhib)>1 | min(inhib)<0) stop("Argument ic must be in interval [0,1].")
  }
  if(normalize!="mean" & normalize!="single"){
    stop(paste(normalize,": Unknown normalization option",sep=""))
  }
  if(graphics!="mean" & graphics!="single" & graphics!="fitted"){
    stop(paste(graphics,": Unknown graphics option",sep=""))
  }
  return(list(measure,control,dilution))
}
