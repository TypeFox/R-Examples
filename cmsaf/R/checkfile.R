checkfile <-
function(infile,outfile){

  pttcheck <- as.logical("FALSE")

  if (exists("ptt")){
    pttcheck <- as.logical("TRUE")
  } else {ptt <- NULL}
    
# check if infile exists, otherwise change path to infile and check again

  status <- as.logical("FALSE") 

  if (!file.exists(infile)&pttcheck){
    file1 <- file.path(ptt,"data",infile)
    if (file.exists(file1)){ 
      outfile <- file.path(ptt,"output",outfile)
      infile <- file1
      status <- as.logical("TRUE")
    }
    if (!file.exists(file1)){
    file1 <- file.path(ptt,"output",infile)
      if (file.exists(file1)){ 
	outfile <- file.path(ptt,"output",outfile)
	infile <- file1
	status <- as.logical("TRUE")
      }else{
	cat(infile," File not found!","\n")
	status <- as.logical("FALSE")
       }
    }
  }else{status <- as.logical("TRUE")}

    status <- status
    infile <- infile
    outfile <- outfile
    output <- list(status,infile,outfile)
    return(output)
}
