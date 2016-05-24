zelig.F <- function(x,fileprefix){
  y <- getwd()
  setwd(x)
  file.copy(paste(fileprefix,"-control-inp.txt",sep=""), "z_control_inp.txt")
  file.copy(paste(fileprefix,"-site-inp.txt",sep=""), "z_site_inp.txt")
  file.copy(paste(fileprefix,"-species-inp.txt",sep=""), "z_species_inp.txt")
  .Fortran("zelig", PACKAGE="seem")
  fileout1 <- paste(fileprefix,"-tracer-out.txt",sep="")
  fileout2 <- paste(fileprefix,"-print-out.txt",sep="")
  file.rename("z_tracer_out.txt", fileout1)
  file.rename("z_print_out.txt", fileout2)
  file.remove("z_control_inp.txt")
  file.remove("z_site_inp.txt")
  file.remove("z_species_inp.txt")
  setwd(y)
  return(paste(x,"/",fileout1,sep="")) 
}

