write.reportProject <- function(mif,mapping){
  if(is.character(mif)){
    data <- read.report(mif,as.list=TRUE)
  } else if (is.list(mif)){
    data <- mif
  }
  else{
    stop("please provide either a path to a mif-file or a mif-file (in list-structure)")
  }
  # read in mapping of the names of variables for the project
  map <- read.csv2(mapping,colClasses="character")
  
  
  
  # select variables and change names of reported variables
  new_data <- list()
  for (n in names(data)){
    for (m in names(data[[n]])){
      new_data[[n]][[m]]<- setNames(data[[n]][[m]][,,map[,names(map)[1]]],map[[names(map)[2]]])
    }
  }
  
  # calculate name of new reporting
  file <- gsub(names(map)[1],names(map)[2],mif)
  # save project reporting
  write.report(new_data,file=file)
}


