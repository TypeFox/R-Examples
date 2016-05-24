#internal function to read WASIM-output files
`myRead.table` <- function (name, subcatchments, has_stat, ts.length, na.values){


   t.data<-array(dim=c(subcatchments+has_stat,ts.length))      #allocate memory for all subbasins and stats column
   col <- paste(sep="", "c", 1:subcatchments)
  if(has_stat==TRUE){
       cols = c("Year","Month","Day","Hour",col, "stat")
  } else {
       cols = c("Year","Month","Day","Hour",col)
  }    

    colClasses=c(rep("integer",4),rep("numeric",subcatchments+has_stat))   #specify data format of column (increases speed, catches "NAN" etc. in output)  


 table <- read.table(name,row.names=NULL, na.strings=na.values, skip=3, nrows=ts.length ,col.names = cols,colClasses=colClasses)
  table[is.na(table)]=NaN     #mark NaNs to be distinguished from NAs        
         nas <- c()
 
 if(NROW(table[[1]])<ts.length){
   cat( "Warning: table ", name, "has length ", NROW(table[[1]]), "expecting length", ts.length, "- filling with NA \n")
   nas <-  rep(NA, ts.length-NROW(table[[1]]))      #NAs to be appended (for filling missing data to achieve desired length)
 }

         for(j in 1:(subcatchments+has_stat)){
          t.data[j,] <- c(table[,j+4], nas)
        }      


    return(t.data)

}
