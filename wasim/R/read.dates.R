#read date information from simulation file
#specify either filename OR nothing

`read.dates` <-
function (file=NULL,
   na.strings=c("999", "999.00","-9999.00"), #strings indicating nodata
   sep="\t",                                  #field delimiter in input file
   skip=3
) {
    
    if (is.null(file))
    {
      filenames=which(data.types$filename!="")
      if (length(filenames)>0)      #use first fully specified filename, if available
      {
        file=data.types$filename[filenames[1]]
      } 
      else
      {
        cat("Error: please use read.dates(file=...) OR specify filename in data.types$filename.\n")
      }
    }
    
    columns<-read.table(file=file,row.names=NULL, skip = skip,nrows=1)  # retrieve columns
    colClasses=rep("NULL",ncol(columns))      #skip all columns...
    colClasses[c(1:4)]="integer"            #...except the first 4
    table <- read.table(file,row.names=NULL, na.strings = na.strings,header=FALSE, skip = skip,colClasses=colClasses)
        
    names(table)= c("Year", "Month", "Day", "Hour")
    return(assemble.date(table))
}

