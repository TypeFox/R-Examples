#Define the write.to.home variable for namespace
write.to.home <- NULL

systemRun <-
  function(directory, output, write.to.home){
    #separate paths for  windows and mac
    if (.Platform$OS.type=="unix"){
      dir <- paste("cd",directory,sep=" ")
      system(paste(dir,"open the_audiolyzr_mac_v5/the_audiolyzr_mac_v5.app", sep="; "))
      
      if (!write.to.home)
        system(paste("cd",output,";","open .."))
    }
    
    if (.Platform$OS.type=="windows"){
      path <- normalizePath(file.path(directory,"the_audiolyzr_win_v5",
                                      "the_audiolyzR_win_v5.exe"))
      shell.exec(path) 
      
      if (!write.to.home)
        shell.exec(file.path(output,".."))
    }
    if (!write.to.home)
      message('To load plots to synthesizer: \n',
              'Drag the folder entitled "json_matrix" \n',
              'to the "Drop folder here to Audiolyze!" box \n',
              'at the bottom of the synthesizer (you may have to scroll down).')
    
    message ("Note: Audio plot files written to \n", file.path(output))
  }

