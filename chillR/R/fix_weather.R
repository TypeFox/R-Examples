#' Weather data fixer and quality checker
#' 
#' This function identifies and interpolates gaps in daily weather records
#' 
#' This function produces a complete record containing all dates between the
#' 1st day of the start year and the last day of the end year (unless the
#' first/last day of the record is after/before these dates - in that case the
#' record is not extended). The values for the columns specified by the columns
#' attribute are linearly interpolated. Missing values during the period
#' indicated by start_date and end_date are added up and summarized in a
#' quality control table.
#' 
#' @param weather a data.frame containing a daily time series dataset. It
#' should have columns c("Year", "Month", "Day") or c("YEAR", "MONTH","DAY") or
#' "YEARMODA".
#' @param start_year integer marking the first year of interest. If not
#' specified, this is assumed to be year 0, which probably means that the
#' entire record will be considered.
#' @param end_year integer marking the last year of interest. If not specified,
#' this is assumed to be year 3000, which probably means that the entire record
#' will be considered.
#' @param start_date start date of the sub-annual period of interest (e.g. the
#' assumed chilling period), defaults to 1 (1st Jan) if not specified
#' @param end_date end date of the sub-annual period of interest (e.g. the
#' assumed chilling period), defaults to 366 (31st Dec, also in non-leap years)
#' if not specified
#' @param columns character vector containing the names of columns of the
#' weather file that should be interpolated and quality checked. If not
#' specified, this defaults to "Tmin" and "Tmax". If these columns don't exist,
#' the function generates an error.
#' @param end_at_present boolean variable indicating whether the interval of
#' interest should end on the present day, rather than extending until the end
#' of the year specified under time_interval[2] (if time_interval[2] is the
#' current year).
#' @return list with two elements: weather: contains the interpolated weather
#' record QC: contains the quality control data.frame, which summarizes missing
#' days, incomplete days (days on which any value is missing), and percentage
#' completeness.
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' fix_weather(KA_weather,2000,2010)
#' 
#' #use a subset of the KA_weather dataset and add an additional day after a gap
#' KA_weather_gap<-rbind(KA_weather,c(Year=2011,Month=3,Day=3,Tmax=26,Tmin=14)) 
#' #fill in the gaps
#' fix_weather(KA_weather_gap, 1990,2011,300,100)
#' 
#' #fix_weather(KA_weather)
#' 
#'  
#' @export fix_weather
fix_weather<-function(weather,start_year=0,end_year=3000,start_date=1,end_date=366,columns=c("Tmin","Tmax"),end_at_present=TRUE)
{
    QC_weather<-function(fixedweather,start_date,end_date,columns,end_at_present)
    {     fixedweather[,"JDay"]<-as.numeric(ISOdate(fixedweather$Year,fixedweather$Month,fixedweather$Day)-ISOdate(fixedweather$Year-1,12,31))
    if(start_date<end_date) {fixedweather[which(fixedweather$JDay>=start_date&fixedweather$JDay<=end_date),"sea"]<-
      fixedweather[which(fixedweather$JDay>=start_date&fixedweather$JDay<=end_date),"Year"]} else
      {fixedweather[which(fixedweather$JDay>=start_date),"sea"]<-
        fixedweather[which(fixedweather$JDay>=start_date),"Year"]+1
      fixedweather[which(fixedweather$JDay<=end_date),"sea"]<-
        fixedweather[which(fixedweather$JDay<=end_date),"Year"]}
    
    if(start_date<end_date) {relevant_days<-start_date:end_date} else
    {relevant_days<-c(start_date:366,1:end_date)}
    seasons<-unique(fixedweather$sea)
    seasons<-seasons[!is.na(seasons)]
 
    if(start_date<1) start_date<-1
    if(end_date<1) end_date<-1
    
    for (sea in seasons)
    {if(leap_year(sea)) {if(end_date>366) end_date_sea<-366 else end_date_sea<-end_date} else 
      {if(end_date>365) end_date_sea<-365 else end_date_sea<-end_date}
     if(leap_year(sea)) {if(start_date>366) start_date_sea<-366 else start_date_sea<-start_date} else
       {if(start_date>365) start_date_sea<-365 else start_date_sea<-start_date}
      if(start_date<end_date) sea_days<-as.numeric(ISOdate(sea-1,12,31)+end_date_sea*86400-(ISOdate(sea-1,12,31)+start_date_sea*86400))+1 else
      sea_days<-as.numeric(ISOdate(sea-1,12,31)+end_date_sea*86400-(ISOdate(sea-2,12,31)+start_date_sea*86400))+1
      
      if (end_at_present) if(sea==as.numeric(format(Sys.Date(),"%Y")))
        if(start_date<end_date) sea_days<-round(as.numeric(min(ISOdate(sea-1,12,31)+end_date_sea*86400,Sys.time())-(ISOdate(sea-1,12,31)+start_date_sea*86400))+1) else
        sea_days<-round(as.numeric(min(ISOdate(sea-1,12,31)+end_date_sea*86400,Sys.time())-(ISOdate(sea-2,12,31)+start_date_sea*86400))+1)

      

    if(sea==seasons[1])
      {gaps<-data.frame(Season=paste(sea-1,"/",sea,sep=""),End_year=sea,Season_days=sea_days,
                        Data_days=length(which(fixedweather$sea==sea)))
       for (ccc in columns) gaps[1,paste("Missing_",ccc,sep="")]<-length(which(fixedweather[which(fixedweather$sea==sea),paste("no_",ccc,sep="")]))+gaps$Season_days[1]-gaps$Data_days[1]
       if(length(columns)==1) gaps[1,"Incomplete_days"]<-max(fixedweather[which(fixedweather$sea==sea),paste("no_",columns,sep="")])+gaps$Season_days[1]-gaps$Data_days[1] else
         gaps[1,"Incomplete_days"]<-max(colSums(fixedweather[which(fixedweather$sea==sea),paste("no_",columns,sep="")]))+gaps$Season_days[1]-gaps$Data_days[1]
       
      # gaps[1,"Incomplete_days"]<-length(which(fixedweather[which(fixedweather$sea==sea),"no_Tmin"]|fixedweather[which(fixedweather$sea==sea),"no_Tmax"]))
                        
      } else
        {df<-data.frame(Season=paste(sea-1,"/",sea,sep=""),End_year=sea,Season_days=sea_days,
                  Data_days=length(which(fixedweather$sea==sea)))
        for (ccc in columns) df[1,paste("Missing_",ccc,sep="")]<-length(which(fixedweather[which(fixedweather$sea==sea),paste("no_",ccc,sep="")]))+df$Season_days[1]-df$Data_days[1]
        if(length(columns)==1) df[1,"Incomplete_days"]<-max(fixedweather[which(fixedweather$sea==sea),paste("no_",columns,sep="")])+df$Season_days[1]-df$Data_days[1] else
          df[1,"Incomplete_days"]<-max(colSums(fixedweather[which(fixedweather$sea==sea),paste("no_",columns,sep="")]))+df$Season_days[1]-df$Data_days[1]
        gaps<-rbind(gaps,df)}
      }
    gaps[,"Perc_complete"]<-round((gaps$Season_days-gaps$Incomplete_days)/gaps$Season_days*100,1)
    
    return(gaps)     
    }      
  
    if((length(names(weather))==2) & ("weather" %in% names(weather))) weather<-weather$weather
    
 fixedweather<-make_all_day_table(weather[which((weather$Year>=start_year)&(weather$Year<=end_year)),])
 if(end_at_present) fixedweather<-fixedweather[which(fixedweather$DATE<Sys.time()),]
 for(ccc in columns)
 {interp<-interpolate_gaps(fixedweather[,ccc])
  fixedweather[,ccc]<-interp$interp
  fixedweather[,paste("no_",ccc,sep="")]<-interp$missing
 }

 return(list(weather=fixedweather,QC=QC_weather(fixedweather,start_date,end_date,columns,end_at_present)))
}
  

  
