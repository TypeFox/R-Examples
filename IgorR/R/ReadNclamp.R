# Functions to summarise the sweep (data) and log files produced by
# Neuromatic/Nclamp

#' Read the log table produced by Nclamp acquisition software for Igor
#' 
#' log tables are special Igor .pxp files that contain only variables.
#' Each entry corresponds to a single run of an Nclamp protocol, 
#' storing information like protocol name, run time etc.  
#' @param f path to the log file.
#' @param Verbose whether to print status information while reading the file.
#' @return A dataframe containing a row for each acquisition protocol run.
#' @author jefferis
#' @export
ReadNclampLogTable<-function(f,Verbose=FALSE){
  x=read.pxp(f,Verbose=Verbose)
  
# WaveFields  
  # Select text fields only
  # cat("orig date field is:",x$vars$FileDate)
  # remove any commas and rationalise spaces
  
  exptDate<-.ParseDate(x$vars$FileDate)
  exptDateStr=format(exptDate)
  if(Verbose) cat("date is:",exptDateStr)
  
  # Do we have a table representation of the data (in waves) ?
  TableFields=sapply(x,function(x) !is.null(attr(x,"WaveHeader")))
  dateFields=NULL
  if(any(TableFields)){
    d=x[TableFields]
    xx=x[TableFields]
  } else{
    selFields=sapply(x[[2]]$vars,class)=="character"
    d=do.call(rbind,lapply(x[-1],function(x) unlist(x$vars[selFields])))
    dateFields=grep("[0-9][0-9]:[0-9][0-9]:[0-9][0-9]",as.character(d[1,]))
    d=data.frame(d)   
  }
  #setdiff(names(x), grep("_[0-9]{3}",names(x),val=T))
  
  for (i in dateFields){
    d[i]=as.POSIXct(paste(exptDateStr,as.character(d[,i])))
#   d[i]=as.POSIXct(strptime(as.character(d[,i]),format="%H:%M:%S"))
  }
  d
}

#' Read all Nclamp log tables from a directory into a list
#' 
#' @param logfiledir path to directory containing log files (pxp files).
#' @param pattern optional regular expression -- see \code{\link{list.files}}.
#' @param ... additional parameters for \code{\link{ReadNclampLogTable}}.
#' @return A named list containing one dataframe for each parsed log file.
#' @author jefferis
#' @export
#' @examples
#' \dontrun{
#' logfiles <- ReadAllNclampLogTables("/GD/projects/PhysiologyData/logs")
#' str(logfiles)
#' }
ReadAllNclampLogTables<-function(logfiledir,pattern="_log[0-9]+[.]pxp$",...){
  logfilenames=dir(logfiledir,pattern=pattern,full.names=T)
  logfiles=list()
  for(i in seq(logfilenames)){
    if(i%%10>0) cat(".") else cat(as.integer(i))
    cat("logfilenames[i]=",logfilenames[i],"\n")
    logfiles[[i]]<-try(ReadNclampLogTable(logfilenames[i]),...)
  }
  names(logfiles)=logfilenames
  logfiles
}

.ParseDate<-function(d)
{
  d=gsub("[, ]+"," ",d)
  exptDate=strptime(d,format="%a %b %d %Y") # "Thu Oct 25 2007"
  if(is.na(exptDate))
    exptDate=strptime(d, format = "%a %e %b %Y") # "Thu 04 Oct 2007"
  if(is.na(exptDate)) stop("Unable to parse date format: ",d)
  exptDate
}

#' Extract summary information from an Nclamp/Igor Sweep File
#' 
#' e.g. for import into Physiology database
#' @param f path to an Nclamp/Igor PXP format sweep file.
#' @param Verbose if \code{TRUE}, print details while parsing underlying PXP
#'   file.
#' @return A list of about 25 fields summarising the sweep file.
#' @author jefferis
#' @export
#' @import tools
#' @examples
#' l <- SummariseSweepFile(system.file("igor", "WedJul407c2_001.pxp", package="IgorR"))
#' cat("There are", l$NumWaves, "waves in the file each of total duration", l$StimWaveLength, 
#'  "ms and sample duration", l$StimSampleInterval, "ms \n") 
SummariseSweepFile<-function(f,Verbose=FALSE){
  s=read.pxp(f,Verbose=Verbose)
  
  fileinfo=file.info(f)
  extraFields=list(FilePath=f,Folder=basename(dirname(f)),
    FileSize=fileinfo$size,FileMTime=fileinfo$mtime,FileMD5=md5sum(f))
  
  chosenFields=c("FileFormat","NumWaves",
    "NumChannels","FileName","FileDate","FileTime","WavePrefix","AcqMode")
  
  summaryFields=s$vars[chosenFields]
  listnames=names(s)[sapply(s,class)=="list"]
  possStimuli=setdiff(listnames, c("vars","Notes",grep("^Chan[A-Z]$",listnames,value=T)))
  if(length(possStimuli)!=1) stop("unable to identify stimulus protocol")
  ProtocolName=possStimuli[1]
  chosenProtocolFields=c("WaveLength","SampleInterval","SamplesPerWave",
    "NumStimWaves","InterStimTime","NumStimReps","InterRepTime","StimRate",
    "RepRate","TotalTime","CurrentFile")
  stimFields=s[[ProtocolName]][["vars"]][chosenProtocolFields]
  names(stimFields)<-paste("Stim",chosenProtocolFields,sep="")
  rval=c(extraFields,summaryFields,ProtocolName=ProtocolName,stimFields)
  nullfields=sapply(rval,is.null)
  rval[nullfields]=NA
  return(rval)
}

#' Summarise multiple sweep files into a single dataframe
#' 
#' Note that this is still a little fragile if the lists produced by 
#' \code{\link{SummariseSweepFile}} do not have consistent field names.
#' @param ff paths to a set of sweep files.
#' @return A dataframe with rows for each sweep file.
#' @author jefferis
#' @seealso SummariseSweepFile
#' @export
SweepFilesToDataFrame<-function(ff){
  ll=lapply(ff,SummariseSweepFile)
  lengths=sapply(ll,length)
  if(any(lengths!=lengths[1]))
    stop("Heterogeneous results from SummariseSweepFile")
  # df=as.data.frame(ll[1],stringsAsFactors=FALSE)
  # if(length(ll)==1) return(df)
  # ll=ll[-1]
  do.call(rbind,lapply(ll,as.data.frame,stringsAsFactors=FALSE))
}

#' Update the CSV file summarising the sweeps in an Nclamp data folder
#' 
#' @param folder path to the folder.
#' @param outfile path to outfile (default: \code{/path/to/datafolder/datafolder.csv}).
#' @param action update newer (default) or force update (not yet implemented).
#' @param DryRun if \code{TRUE}, report which files would be processed, but do
#'   nothing.
#' @return \code{TRUE} if something happened, \code{FALSE} otherwise.
#' @author jefferis
#' @export
UpdateSweepDataFrame<-function(folder,outfile=NULL,action=c("update","force"),DryRun=FALSE){
  action=match.arg(action)
  folder=path.expand(folder) # replace ~ by full path if necessary
  if(is.null(outfile)){
    outfile=file.path(folder,paste(basename(folder),sep=".","csv"))
  }
  infiles=dir(folder,pattern="_[0-9]+\\.pxp$",full.names=T)
  if(!file.exists(outfile)){
    # we need to process all infiles
    newinfiles=infiles
  } else {
    intimes=file.info(infiles)$mtime
    outtime=file.info(outfile)$mtime
    newinfiles=infiles[intimes>outtime]
  }
  if(length(newinfiles)>0){
    if(DryRun) {
      cat("Would process the following files:\n",sep="",paste(newinfiles,collapse="\n"))
    } else {
      newdf=SweepFilesToDataFrame(newinfiles)
      if(file.exists(outfile)){
        # read in the old data frame, and overwrite updated rows/append new rows
        olddf=read.csv(outfile,stringsAsFactors=FALSE)
        rownames(olddf)=olddf[,1] # set the rownames to first col which should be FilePath
        olddf[rownames(newdf),]=newdf
        newdf=olddf
      }
      write.csv(newdf,outfile)
      return(invisible(TRUE))
    }
  }
  return(invisible(FALSE))
}
