ncdump <- function(filename, args = "-h"){
  cmd <- paste("ncdump", args, filename, sep = " ")

  if(is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
    if(pbdMPI::comm.rank() == 0){
      system(cmd)
    }
  } else{
    system(cmd)
  }

  invisible()
}
