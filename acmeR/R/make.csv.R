#'@import utils
#'
#'
make.csv <- function(dframe, fpath, fname){
  if(missing(dframe)){
    stop("Missing data frame to convert");
  }else if(!(class(dframe) %in% c("data.frame","matrix"))){
    stop("Requires a variable of type data.frame or matrix")
  }
  if(missing(fpath)){
    fpath <- getwd();
  }
  if(missing(fname)){
    fname <- deparse(substitute(dframe))
  }
  len <- nchar(fname)
  if(substr(fname, len-4+1, len) != ".csv"){
    fname <- paste(fname,".csv",sep="")
  }
  
  write.csv(dframe, file = paste(fpath,"/",fname,sep=""));
}