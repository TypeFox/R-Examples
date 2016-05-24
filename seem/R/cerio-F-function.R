cerio.F <- function(x,fileprefix){
  y <- getwd()
  setwd(x)
  file.copy(paste(fileprefix,"-inp.txt",sep=""), "cerio_inp.txt")
  .Fortran("cerio", PACKAGE="seem")
  fileout <- paste(fileprefix,"-out.txt",sep="")
  filedex <- paste(fileprefix,"-dex.txt",sep="")
  file.rename("cerio_out.txt", fileout)
  file.rename("cerio_dex.txt", filedex)
  file.remove("cerio_inp.txt")
  setwd(y)
  return(paste(x,"/",fileout,sep="")) 
}

