forsucc.C <- function(x,fileprefix){
  y <- getwd()
  setwd(x)
  file.copy(paste(fileprefix,"-inp.txt",sep=""), "forsucc_inp.txt")
  .C("forsucc", PACKAGE="seem")
  fileout <- paste(fileprefix,"-out.txt",sep="")
  filedex <- paste(fileprefix,"-dex.txt",sep="")
  file.rename("forsucc_out.txt", fileout)
  file.rename("forsucc_dex.txt", filedex)
  file.remove("forsucc_inp.txt")
  setwd(y)
  return(paste(x,"/",fileout,sep="")) 
}

