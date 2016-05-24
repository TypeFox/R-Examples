dataPreproc<-function(dataDir=getwd(), blocks=12, spot="aushon", exportNo=3, correct="both", remove_flagged=NULL){
                   
  ################################################################################
  # 1. Import and convert raw data from ".gpr"-files, slide- & sampledescription #
  ################################################################################
  
  setwd(dataDir) 
  rawdat<-read.Data(blocksperarray=blocks, spotter=spot, remove_flagged=remove_flagged) # list of length 4, see read.Data()
  
  # create an analysis folder labeled by the date of analysis
  DIR <- paste("analysis",sub(" .*$","",Sys.time()),sep="_")
  if(!file.exists(DIR)){
    dir.create(DIR)
  }
  
  # write data to Excel sheet format, see write.Data()
  setwd(DIR)
  write.Data(rawdat)
  # remove NA-rows in RPPanalyzer outputs (i.e. no measurement)
  fgRaw.tmp<-read.delim("Dataexpression.txt", stringsAsFactors=FALSE, row.names=NULL, header=TRUE) # foreground intensities
  fgRaw<-read.delim("Dataexpression.txt", 
                    skip=max(which(fgRaw.tmp[,1]==""))+1, stringsAsFactors=FALSE, row.names=NULL, header=TRUE)
  bgRaw.tmp<-read.delim("Databackground.txt", stringsAsFactors=FALSE, row.names=NULL, header=TRUE) # background intensities
  bgRaw<-read.delim("Databackground.txt", 
                    skip=max(which(bgRaw.tmp[,1]==""))+1, stringsAsFactors=FALSE, row.names=NULL, header=TRUE)
  fgNAVec<-which(is.na(fgRaw[,"ID"]))
  bgNAVec<-which(is.na(bgRaw[,"ID"]))
  if(length(fgNAVec)>0){
    fgRaw<-fgRaw[-fgNAVec,] 
  }
  if(length(bgNAVec)>0){
    bgRaw<-bgRaw[-bgNAVec,]
  }
  colnames(fgRaw)<-sub("X","",gsub("\\.","-",colnames(fgRaw)))
  colnames(bgRaw)<-sub("X","",gsub("\\.","-",colnames(bgRaw)))
  
  arrayDesc<-rawdat$arraydescription
  
  
  ############################################################
  # 2. FOLD CHANGE calculation according to correctDilinterc #
  ############################################################
  
  if(correct!="none"){ # correct intercepts for targets (and FCF)
    
    # split raw data (fgRaw) into dilution data (dilData) and measurement data to be normalized
    dilData <- fgRaw[which(fgRaw$sample_type=="control" & !is.na(fgRaw$dilSeriesID)),]
    
    # warning in case of exportNo = 4 and only one dilseriesID value
    if(length(unique(dilData$dilSeriesID))==1 & exportNo == 4)
    {cat("warning: exportNo has to be reduced or different dilseriesID flags have to be used.")}
    
    # adapt signal intensities via subtraction of dilution intercept at concentration 0
    normdatFC<-correctDilinterc(dilseries=dilData, arraydesc=arrayDesc, timeseries=fgRaw[which(fgRaw$sample_type=="measurement"),], 
                                exportNo=exportNo)
    
    # correct data values for negative ones
    filedata = c()
    for(A1 in colnames(arrayDesc))  
    { ind = arrayDesc["array.id",which(colnames(arrayDesc) ==A1)]
      if (min(normdatFC[, ind]) < 0) 
      { filedata = rbind(filedata, c(A1, arrayDesc["target",A1], abs(min(normdatFC[, ind]))) )
        normdatFC[, ind] <- normdatFC[, ind] + abs(min(normdatFC[, ind])) +1
        colnames(filedata) = c("array ID", "target", "(slide+pad)-specific abs(min)")}
    }
    if(!is.null(filedata))
    {write.table(filedata, file = "negval_correction.txt", sep = "\t", row.names = F)
    }else{
    cat("no negative value correction was necessary\n")}
        
    # get FCF columns
    fcfCol<-colnames(arrayDesc)[grep("protein",arrayDesc["target",])]
    
    # print warning, if FCF intercept > 50% FCF signal
    fcfInterc<-fgRaw[which(fgRaw$sample_type=="measurement"),fcfCol]-normdatFC[,fcfCol]
    if(any(fcfInterc > 0.5*fgRaw[which(fgRaw$sample_type=="measurement"),fcfCol])){
      cat("warning: FCF correction factor exceeds 50% of FCF signal\n")
    }
    
    if(correct=="noFCF"){ # reset old FCF values
      normdatFC[,fcfCol]<-fgRaw[which(fgRaw$sample_type=="measurement"),fcfCol]
    }
    cordat<-list()
    cordat[[1]]<-as.matrix(normdatFC[,colnames(arrayDesc)])
    #dummy matrix, BG neglectible --> not corrected for intercepts
    cordat[[2]]<-as.matrix(bgRaw[which(bgRaw$sample_type=="measurement"),colnames(cordat[[1]])])
    cordat[[3]]<-rawdat$arraydescription[,colnames(cordat[[1]])]
    cordat[[4]]<-rawdat$sampledescription
    if(length(fgNAVec)>0){
      cordat[[4]]<-rawdat$sampledescription[-fgNAVec,]
    }
    cordat[[4]]<-cordat[[4]][which(cordat[[4]]$sample_type=="measurement"),]
    names(cordat)<-names(rawdat)
  }else{ # do not use intercept correction at all
    cordat <- list()
    cordat[[1]] <- as.matrix(fgRaw[which(fgRaw$sample_type == "measurement"), colnames(arrayDesc)])
    cordat[[2]] <- as.matrix(bgRaw[which(bgRaw$sample_type == "measurement"), colnames(cordat[[1]])])
    cordat[[3]] <- rawdat$arraydescription[, colnames(cordat[[1]])]
    cordat[[4]] <- rawdat$sampledescription
    if (length(fgNAVec) > 0) {
      cordat[[4]] <- rawdat$sampledescription[-fgNAVec, ]
    }
    cordat[[4]] <- cordat[[4]][which(cordat[[4]]$sample_type == "measurement"), ]
    names(cordat) <- names(rawdat)
  }
  
  
  ########################
  # 3. FCF normalization #
  ########################
  
  # normalizeRPPA() with method "proteinDye" uses FCF signal for normalization 
  # the median of the normalizer values is added after normalization of log2 data and the data is returned at native scale
  normdat<-normalizeRPPA(cordat, method="proteinDye", vals="native") # list of length 4, see read.Data()
    
  
  ##############################
  ## 4. Quality Control plots ##
  ##############################
  
  # plot target and blank signal from serially diluted control samples of raw RPPA data set
  plotQC(rawdat,file="QC_dilutioncurve_raw.pdf")
  
  # plot the blank signals and the target specific signals of dilution intercept corrected and FCF normalized RPPA data
  plotMeasurementsQC(normdat, file = "QC_targetVSblank_normed.pdf", arrays2rm = c("protein"))
 
  # qq-plot
  plotqq(normdat, fileName = "QC_qqPlot_normed.pdf")
  
  setwd(dataDir)
  return(list(rawdat=rawdat, cordat=cordat, normdat=normdat, DIR=DIR))
}