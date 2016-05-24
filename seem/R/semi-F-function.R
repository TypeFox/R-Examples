semi.F <- function(x,fileprefix){
  y <- getwd()
  setwd(x)
  file.copy(paste(fileprefix,"-inp.txt",sep=""), "semi_inp.txt")
  .Fortran("semi", PACKAGE="seem")
  fileout <- paste(fileprefix,"-out.txt",sep="")
  file.rename("semi_out.txt", fileout)
  file.remove("semi_inp.txt")
  setwd(y)
  return(paste(x,"/",fileout,sep="")) 
}

