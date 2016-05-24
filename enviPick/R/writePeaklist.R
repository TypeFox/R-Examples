writePeaklist <-
function(MSlist,directory,filename,overwrite=FALSE){

      ##########################################################################
      if(!length(MSlist)==8){stop("This is not an MSlist object")}
      if(!MSlist[[1]][[5]]){stop("MSlist does not contain picked peaks - abort.")}
      if(!file.exists(directory)){stop("invalid directory")}     
      fileout<-paste(directory,"\\",filename,sep="")
      if(file.exists(fileout) & !overwrite){stop("file already exists; cannot overwrite!")}
      ##########################################################################      
      write.table(MSlist[[8]],file=fileout)
      ##########################################################################

}
