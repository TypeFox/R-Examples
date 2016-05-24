saveMarkerModels <-
function(
            markers=NA, data, diplo=NA,
            select=TRUE, diploselect=TRUE, # by default all samples of tetraploids and diploids selected
            maxiter=40, maxn.bin=200, nbin=200,
            sd.threshold=0.1, p.threshold=0.99,
            call.threshold=0.6, peak.threshold=0.85,
            try.HW=TRUE, dip.filter=1,
            sd.target=NA,
            ncores=NA,
            logfile="", modelfile, allmodelsfile="", 
            scorefile, diploscorefile="",
            plot="none", plot.type="png") {
  if (class(data)!="data.frame" || 
      nrow(data)==0 ||  
      length(which(names(data)=="MarkerName"))!=1 ||
      length(which(names(data)=="SampleName"))!=1 ||
      length(which(names(data)=="ratio"))!=1) {
    stop("data is not a valid data frame")
  } else {
    if (length(select)!=nrow(data)) {
      select <- rep(select, times=ceiling(nrow(data)/length(select)))
      if (length(select)>nrow(data)) select <- select[1:nrow(data)]
    }
  } 
  if (class(diplo)=="data.frame") {
    if (nrow(diplo)==0 ||
        length(which(names(diplo)=="MarkerName"))!=1 ||
        length(which(names(diplo)=="SampleName"))!=1 ||
        length(which(names(diplo)=="ratio"))!=1) {
      stop("diplo is not a valid data frame")
    } else {
      if (length(diploselect)!=nrow(diplo)) {
        diploselect <- rep(diploselect, times=ceiling(nrow(diplo)/length(diploselect)))
        if (length(diploselect)>nrow(diplo)) diploselect <- diploselect[1:nrow(diplo)]
      }
    }  
  }
  plot <- tolower(plot)
  if (plot!="none") {
    plot.type <- check.plottype(plot.type)
    if (plot.type=="none") plot<-"none"  
  }
  if (plot=="none") plot.dir <- NA else plot.dir <- getPlotsDirectory(logfile)
  plot.type <- paste("_",plot.type,sep="")
  markernames <- levels(as.factor(as.character(data$MarkerName)))
  samplenames <- levels(as.factor(as.character(data$SampleName)))
  diplonames <- character(0)
  if (class(diplo)=="data.frame") { #not NA
    diplonames <- levels(as.factor(as.character(diplo$SampleName)))
  }  
  if (length(diplonames)==0) diploscorefile<-""
  if (length(markers)==1 && is.na(markers)) markers <- 1:length(markernames)
  dip.filter <- as.integer(dip.filter) # for compatibility with earlier versions where it was logical
  suppressWarnings( {
    file.remove(scorefile)
    file.remove(diploscorefile)
    file.remove(allmodelsfile)
    file.remove(modelfile)
    file.remove(logfile)
    }
  )
  if (logfile!="") {
    write (paste("Log-file started at",format(Sys.time(),format="%Y%m%d-%H%M%S")),file=logfile)
    write (paste("sample count in data =",length(samplenames)),file=logfile,append=T)
    write (paste("sample count in diplo =",length(diplonames)),file=logfile,append=T)
    write (paste("marker count in data =",length(markernames)),file=logfile,append=T)
    write (paste("marker count to score =",length(markers)),file=logfile,append=T)
    write (paste("maxiter =",maxiter),file=logfile,append=T)
    write (paste("maxn.bin =",maxn.bin),file=logfile,append=T)
    write (paste("nbin =",nbin),file=logfile,append=T)
    write (paste("sd.threshold =",sd.threshold),file=logfile,append=T)
    write (paste("p.threshold =",p.threshold),file=logfile,append=T)
    write (paste("call.threshold =",call.threshold),file=logfile,append=T)
    write (paste("peak.threshold =",peak.threshold),file=logfile,append=T)
    write (paste("try.HW =",try.HW),file=logfile,append=T)
    write (paste("dip.filter =",dip.filter),file=logfile,append=T)
    write (paste("sd.target =",sd.target),file=logfile,append=T)
  }
  if (length(ncores)!=1 || is.na(ncores) || ncores<2 || .Platform$OS.type!="unix") ncores <- 1
  if (ncores>1 && !(require("doMC",quietly=T) && require("foreach",quietly=T))) ncores <- 1
  if (ncores==1) {
    for (mrk in markers) if (mrk %in% 1:length(markernames)) {
      mrkresult <- fitTetra(
              mrk, data, diplo,
              select, diploselect,
              maxiter, maxn.bin, nbin,
              sd.threshold,
              p.threshold,
              call.threshold,
              peak.threshold,
              try.HW,
              dip.filter,
              sd.target,
              plot, plot.type, plot.dir) 
      write.files(mrkresult, modelfile, scorefile, diploscorefile, logfile, allmodelsfile)           
    } # for mrk
  } else { #ncores>1  
    registerDoMC(cores=ncores)
    writing.files <- F #to avoid collisions between threads
    foreach (mrk = markers) %dopar% {
      mrkresult <- fitTetra(
              mrk, data, diplo,
              select, diploselect,
              maxiter, maxn.bin, nbin,
              sd.threshold,
              p.threshold,
              call.threshold,
              peak.threshold,
              try.HW,
              dip.filter,
              sd.target,
              plot, plot.type, plot.dir)
      while (writing.files) {
        Sys.sleep(0.5)
      }  
      writing.files <- T 
      write.files(mrkresult, modelfile, scorefile, diploscorefile, logfile, allmodelsfile)
      writing.files <- F     
    } # foreach mrk
    #sort output files:
    if (file.exists(modelfile)) {
      dat <- read.table(modelfile,header=TRUE,na.strings="",sep="\t")
      dat <- dat[order(dat$marker),]
      write.table(dat,file=modelfile,sep="\t",na="",row.names=F,col.names=T,quote=F)
    }
    if (file.exists(scorefile)) {  
      dat <- read.table(scorefile,header=TRUE,na.strings="",sep="\t")
      dat <- dat[order(dat$marker,dat$sample),]
      write.table(dat,file=scorefile,sep="\t",na="",row.names=F,col.names=T,quote=F)
    }
    if (allmodelsfile!="" && file.exists(allmodelsfile)) {  
      dat <- read.table(allmodelsfile,header=TRUE,na.strings="",sep="\t")
      dat <- dat[order(dat$marker,dat$m),]
      write.table(dat,file=allmodelsfile,sep="\t",na="",row.names=F,col.names=T,quote=F)
    }  
  }
  Sys.sleep(5) #allow a few sec to save all files in batch mode
}
