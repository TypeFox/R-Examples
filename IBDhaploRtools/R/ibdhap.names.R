## function to return names of chromosomes in each set

ibdhap.names <- function( ids.file = NULL, ids.filename = NULL){

  ## information from ids file
  if( !is.null(ids.file)){ ## specified the ids file already loaded

    par.file <- ids.file
    print("ids matrix accepted")
    
  }else if(is.null(ids.file)){ ## specified a file location to read
    
    n.col.par <- max(count.fields(ids.filename)) ## max number of columns in the input file
    par.file <- read.table(ids.filename, fill =TRUE, colClasses = "numeric", col.names=1:n.col.par)
    print("ids matrix read from file")

  }else{ stop("Invalid ids input") }

  name.cols <- seq(4, ncol(par.file)) ## which columns are identifiers
  
  name.list <- par.file[, name.cols] ## collect identifiers
  rownames(name.list) <- par.file[,1] ## change row name to set name
  
  return(name.list)
}


  
