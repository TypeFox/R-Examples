read.multiple.objects <- function (directory, prefix, class=c("spc", "vgc", "tfl")) {

  ## this is the R way of checking an option against a fixed range of values
  class <- match.arg(class)

  # obtain list of files matching pattern in specified directory and
  # complain and stop if there is no file matching pattern
  files <- list.files(directory,pattern=paste("^",prefix,"[\\.]",".+","[\\.]",class,"$",sep=""))
  file.number = length(files)
  if (file.number==0) {
    stop("no file matching pattern ",prefix,".*.",class," found")
  }

  # now, we go through the list of files and read them to a list
  index <- 1
  results<-list()
  
  while (index <= file.number) {

    # extract whatever occurs between the common prefix and the
    # common suffix, using it as an id for the current object
    id <- sub(paste("^",prefix,"\\.",sep=""),"",files[index])
    id <- sub(paste("\\.",class,"$",sep=""),"",id)

    # depending on object type, use different read functions
    if (class == "spc") {
      results[[id]] <- read.spc(paste(directory,files[index],sep="/"))
    }
    else {
      if (class == "vgc") {
        results[[id]] <- read.vgc(paste(directory,files[index],sep="/"))
      }
      else{
        results[[id]] <- read.tfl(paste(directory,files[index],sep="/"))
      }
    }
    
    index <- index + 1
      
  }

  return(results)
  
}
