audiolyzRcheck <- function(directory){
  #Check for audiolyzR installation
  if (.Platform$OS.type=="windows"){
    if(!file.exists(paste(directory,"the_audiolyzr_win_v5",sep="/"))) {
      ANSWER <- readline("The latest audiolyzR standalone synthesizer does not appear to be installed. \nWould you like to download and install the latest Windows version now? \nType y or n")
      if (substr(ANSWER,1,1)=="y"){
        if(!file.exists(paste(directory,"the_audiolyzr_win_v5",sep="/")))
          dir.create(paste(directory,"the_audiolyzr_win_v5",sep="/"), recursive=TRUE, showWarnings=FALSE)
        tmpdir <- tempdir()
        download.file("http://s3.amazonaws.com/audiolyzR/installers/the_audiolyzR_win_v5.zip",
                      destfile=paste(tmpdir, "the_audiolyzr_win_v5.zip", sep="/"))
        unzip(paste (tmpdir,"the_audiolyzr_win_v5.zip",sep="/"),
              exdir = directory)
        cat("NOTE: The synthesizer has been installed in",
            paste(paste(directory,"the_audiolyzr_win_v5",sep="/")))
    } else {
      stop("Please see ?audiolyzR for details about how to install the external synthesizer. \nOnce this is completed, please try your function again.")
    }
    }
    #kills an already-running instance of the synthesizer
    options(warn=-1)
    shell("taskkill /IM the_audiolyzr_win_v5.exe /F >nul 2>nul", ignore.stdout=TRUE, ignore.stderr=TRUE)
    options(warn=0)
  }
     
  if (.Platform$OS.type=="unix"){
    #check whether the audio app is installed
    if(!file.exists(paste(directory,"/the_audiolyzr_mac_v5/the_audiolyzr_mac_v5.app",sep=""))) 
    {
      ANSWER <- readline("The latest audiolyzR synthesizer does not appear to be installed. \nWould you like to download the latest Mac version now? \n type y or n")
      if (substr(ANSWER,1,1)=="y"){
        if(!file.exists(paste (directory,"the_audiolyzr_mac_v5",sep="/")))
          dir.create(paste (directory,"the_audiolyzr_mac_v5",sep="/"), recursive=TRUE, showWarnings=FALSE)
        download.file("http://s3.amazonaws.com/audiolyzR/installers/the_audiolyzR_mac_v5.zip",
                      destfile=paste(directory, "the_audiolyzr_mac_v5", "the_audiolyzr_mac_v5.zip", sep="/"))
        system(paste("unzip", paste(directory, "the_audiolyzr_mac_v5", "the_audiolyzr_mac_v5.zip", sep="/"),
                     "-d", paste (directory,"the_audiolyzr_mac_v5",sep="/"), sep=" "), 
               ignore.stdout=TRUE)
        message("NOTE: The synthesizer has been installed in ",
                paste(directory,"the_audiolyzr_mac_v5",sep="/")) 
      } else {
        stop("Please see ?audiolyzR for details about the external synthesizer. \nOnce installation is completed, please try your function again.")
      }
    }
    system("killall -9 the_audiolyzR_mac_v5", ignore.stderr=TRUE)
  }
}