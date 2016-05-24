PWbatch <-
function(folderin,folderout,msconvert_path,notintern=FALSE,use_format="mzXML"){

      ##########################################################################
      # checks & setups ########################################################
      if(nchar(Sys.which("msconvert")[[1]])==0){
        cat("msconvert not in system path - ok if msconvert_path correct")
      }
      if(substr(folderin,nchar(folderin)-1,nchar(folderin))!="\\"){
        folderin<-paste(folderin,"\\",sep="")
      }
      that<-shell(shQuote(msconvert_path),intern=TRUE)
      if(length(that)!=149){
        stop(" -> but: invalid msconvert_path!")
      }
      files<-list.files(folderout);
      files<-files[files!="temp"];
      if(length(files)>0){
        stop("folderout not empty! Define empty folder!")
       }
      ##########################################################################
      # convert to mzML ########################################################
      files<-list.files(folderin);
      prog<-winProgressBar("Convert .Raw to .mzXML...",min=0,max=length(files));
      setWinProgressBar(prog, 0, title = "Convert .Raw to .mzXML...", label = NULL)
      for(i in 1:length(files)){
        setWinProgressBar(prog, i, title = "Convert .Raw to .mzXML...", label = NULL)
        if(
          substr(files[i],nchar(files[i])-3,nchar(files[i]))==".RAW" ||
          substr(files[i],nchar(files[i])-3,nchar(files[i]))==".raw"
        ){
          there1<-paste(folderin,files[i],sep="");
          there2<-paste(" -o ",shQuote(folderout),sep="");
		  filtered0<-paste(shQuote("--"),use_format,sep="")
          filtered1<-paste(" --filter ",shQuote("peakPicking true 1"),sep="")
          filtered2<-paste(" --filter ",shQuote("msLevel 1"),sep="")
          system(
              paste(
                shQuote(msconvert_path),
                shQuote(there1),
				filtered0,
                filtered1,
                filtered2,
                there2
              )
          ,intern=notintern)
        }else{
          cat("Warning: Found not .RAW file in folderin - file skipped!\n")
        }
      }
      close(prog);
      ##########################################################################

}
