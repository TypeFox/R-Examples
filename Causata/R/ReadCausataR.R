require("R.utils") # used for the countLines function

ReadCausataR <- function(rFile, countRows=FALSE){
  # reads an R file exported with Causata data and returns the variable names and types
  # extract column classes
  scanStr <- scan(file=rFile, what=character(), sep="\n", quiet=TRUE)
  
  # the input file should have 1 or 2 rows, if 2 then focus there and disregard first
  if (length(scanStr) == 2){
    scanStr <- scanStr[[2]] # use the second row
  }
  split1 <- strsplit(scanStr, "colClasses=c(", fixed=TRUE)
  split2 <- strsplit(split1[[1]][2], ")", fixed=TRUE)
  colClasses <- eval( parse(text=paste("c(", split2[[1]][1], ")")) )
  
  # extract column names
  split1 <- strsplit(scanStr, "col.names=c(", fixed=TRUE)
  split2 <- strsplit(split1[[1]][2], ")", fixed=TRUE)
  col.names <- eval( parse(text=paste("c(", split2[[1]][1], ")")) )
  
  # extract CSV file name
  # the filename as exported from causata does not include the path
  split1 <- strsplit(scanStr, "read.csv(", fixed=TRUE)
  split2 <- strsplit(split1[[1]][2], ",", fixed=TRUE)
  fileData <- eval( parse(text=split2[[1]][1]) )
  # add on the path provided in rFile
  fileData <- file.path(dirname(rFile), fileData)
  
  # ensure that the variable names and column classes have the same length
  stopifnot(length(colClasses) == length(col.names))
  
  # count the number of lines in the data file
  # this uses a function from the R.utils package
  # this can take a while for larger files
  if (countRows){
    nrows <- countLines(fileData)
  } else {
    nrows <- -1 # initialize to -1, negative values are ignored by read.csv  
  }
  
  # build a list of outputs
  outlist <- list(
    fileR = rFile, # filename of the R script
    fileData = fileData, # filename of the csv data
    colClasses = colClasses, # column classes
    col.names = col.names, # column names
    nrows = nrows
  )
  return(outlist)
}
