river.C <- function(x,fileprefix){
  y <- getwd()
  setwd(x)
  file.copy(paste(fileprefix,"-inp.txt",sep=""), "river_inp.txt")
  .C("river", PACKAGE="seem")
  fileout <- paste(fileprefix,"-out.txt",sep="")
  filedex <- paste(fileprefix,"-dex.txt",sep="")
  file.rename("river_out.txt", fileout)
  file.rename("river_dex.txt", filedex)
  file.remove("river_inp.txt")
  setwd(y)
  return(paste(x,"/",fileout,sep="")) 
}

