#' Fill in missing days in incomplete time series
#' 
#' Time series often have gaps, and these are often not marked by 'no data'
#' values but simply missing from the dataset. This function completes the time
#' series by adding lines for all these missing records. For these lines, all
#' values are set to 'NA'.
#' 
#' 
#' @param tab a data.frame containing a time series dataset. It should have
#' columns c("Year", "Month", "Day") or c("YEAR", "MONTH","DAY") or "YEARMODA".
#' @return data frame containing all the columns of the input data frame, but
#' one row for each day between the start and end of the dataset. Data values
#' for the missing rows are filled in as 'NA'. Dates are expressed as
#' c("YEARMODA","DATE","Year","Month","Day"). In this, 'DATE' is the date in
#' ISOdate format.
#' @author Eike Luedeling
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' @keywords utility
#' @examples
#' 
#' #use a subset of the KA_weather dataset and add an additional day after a gap
#' KA_weather_gap<-rbind(KA_weather[1:100,],c(Year=1998,Month=6,Day=3,Tmax=26,Tmin=14)) 
#' #fill in the gaps
#' make_all_day_table(KA_weather_gap)
#' 
#' 
#' @export make_all_day_table
make_all_day_table<-function(tab) #tab should have columns named Year, Month and Day (or YEAR, MONTH, DAY; or YEARMODA)
{
 columns<-colnames(tab)  
  if("YEAR" %in% colnames(tab)) colnames(tab)[which(colnames(tab)=="YEAR")]<-"Year"
  if("year" %in% colnames(tab)) colnames(tab)[which(colnames(tab)=="year")]<-"Year"
  if("MONTH" %in% colnames(tab)) colnames(tab)[which(colnames(tab)=="MONTH")]<-"Month"
  if("month" %in% colnames(tab)) colnames(tab)[which(colnames(tab)=="month")]<-"Month"
  if("DAY" %in% colnames(tab)) colnames(tab)[which(colnames(tab)=="DAY")]<-"Day"
  if("day" %in% colnames(tab)) colnames(tab)[which(colnames(tab)=="day")]<-"Day"
  if("YEARMODA" %in% colnames(tab))
      {tab[,"Year"]<-trunc(tab[,"YEARMODA"]/10000)
       tab[,"Month"]<-trunc((tab[,"YEARMODA"]-tab[,"Year"]*10000)/100)
       tab[,"Day"]<-tab[,"YEARMODA"]-tab[,"Year"]*10000-tab[,"Month"]*100
      }
tab[,"YEARMODA"]<-as.numeric(tab[,"Year"])*10000+as.numeric(tab[,"Month"])*100+as.numeric(tab[,"Day"])
tab[,"DATE"]<-ISOdate(tab[,"Year"], tab[,"Month"], tab[,"Day"])
datevec<-seq(tab$DATE[1], tab$DATE[nrow(tab)], "day")
tab2<-tab[,!(names(tab) %in% c("Year","Day","Month","Day","DATE"))]

output<-data.frame(DATE=datevec)
output[,"Year"]<-as.numeric(format(output[,"DATE"], "%Y"))
output[,"Month"]<-as.numeric(format(output[,"DATE"], "%m"))
output[,"Day"]<-as.numeric(format(output[,"DATE"], "%d"))
output[,"YEARMODA"]<-output[,"Year"]*10000+output[,"Month"]*100+output[,"Day"]

output<-merge(output,tab2,by="YEARMODA",all.x=TRUE)
#output<-output[,c(2,3,4,5,11,12)]
#colnames(output)<-c("DATE","Year","Month","Day","Tmin","Tmax")
return(output)
}

#KA_weather_mod<-KA_weather[1:100,]
#KA_weather_mod<-rbind(KA_weather_mod,c(Year=1998,Month=6,Day=3,Tmax=26,Tmin=14))
#tab<-KA_weather_mod
#make_all_day_table(KA_weather_mod)
