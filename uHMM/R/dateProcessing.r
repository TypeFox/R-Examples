# Fonction dateProcessing
#' @title Create a vector of observation times
#' @description Paste date and time of observations in a vector.
#' @param data dataframe with a column "Dates" in "YYYY-MM-DD" format and a column "Hours" in "HH:MM:SS" format.
#' @return numeric vector (use \code{\link[chron]{chron}} function to get the vector in date format).
#' @importFrom chron chron


.dateProcessing<-function(data){
  
  year<-substring(data[,"Dates"],1,4)
  month<-substring(data[,"Dates"],6,7)
  day<-substring(data[,"Dates"],9,10)
  time<-data[,"Hours"]
  seconds<-as.data.frame(table(substring(time,7,8)))
  if (nrow(seconds)==1 & seconds[1,1]=="00"){
    time<-paste(substr(time,1,6),"01",sep="")
  }
  
    return(as.numeric(chron(
        paste(day,month,year,sep="/"),time,format=c("d/m/y","h:m:s")
        )))
}