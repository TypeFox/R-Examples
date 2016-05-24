getPlotsDirectory <-
function(logfile) {
  if (!is.na(logfile)) 
    logfile <- gsub("(^ +)|( +$)", "", logfile) #strip leading and trailing blanks
  if (is.na(logfile) || logfile=="") plotdir <- "plots"
  else {
    #if logfile ends in ".log" remove the extension:
    i<-nchar(logfile)
    while (i>0 && substring(logfile,i,i)!=".") i<-i-1;
    if (i==0) plotdir <- paste("plots_",logfile,sep="")
    else {
      if (tolower(substring(logfile,i,nchar(logfile)))==".log") {
        logfile <- substring(logfile,1,i-1)
      }  
      if (logfile=="") plotdir <- "plots" 
      else plotdir <- plotdir <- paste("plots_",logfile,sep="")
    }
  }  
  dir.create(plotdir)
  #return value:
  if (file.exists(plotdir)) paste(plotdir,"/",sep="")
  else ""
}
