PWfile <-
function(rawfile,folderout,msconvert_path,notintern=FALSE,use_format="mzXML"){

      ##########################################################################
      # checks & setups ########################################################
      if(nchar(Sys.which("msconvert")[[1]])==0){
        cat("msconvert not in system path - ok if msconvert_path correct")
      }
      if(
          sum(substr(rawfile,nchar(rawfile)-3,nchar(rawfile))!=".RAW",substr(rawfile,nchar(rawfile)-3,nchar(rawfile))!=".raw")!=1
      ){stop("rawfile not a .RAW file")}	  
      ##########################################################################
      # convert ################################################################
      there2<-paste(" -o ",shQuote(folderout),sep="");
	  filtered0<-paste(shQuote("--"),use_format,sep="")
      filtered1<-paste(" --filter ",shQuote("peakPicking true 1"),sep="")
      filtered2<-paste(" --filter ",shQuote("msLevel 1"),sep="")
      system(
              paste(
                shQuote(msconvert_path),
                shQuote(rawfile),
				filtered0,
                filtered1,
                filtered2,
                there2
              )
      ,intern=notintern)
      ##########################################################################

}
