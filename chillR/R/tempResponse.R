#' Calculation of climatic metrics from hourly temperature records
#' 
#' Extension of the chilling function, which calculated four pre-defined
#' temperature-based metrics. This function has more flexibility, because it
#' allows specifying the models that should be calculated. These can be
#' selected from a small set of models provided with chillR, but they can also
#' be defined by the user. Precondition at the moment is that they require
#' hourly temperature only as inputs.
#' 
#' The function calculates the total of user-specified temperature-based
#' metrics over periods delineated by Start_JDay and End_JDay. Models for
#' calculating these metrics are provided in the models list, whose elements
#' are named functions that convert hourly temperature records into a
#' cumulative record of the climate metric of interest. The metric is then
#' added up cumulatively over the entire temperature record and then summarized
#' by season. Examples of functions that can be used are Chilling_Hours,
#' Utah_Model, Dynamic_Model and GDH. The custom_model function allows
#' customized simply weight-based models, which assign differential weights to
#' temperatures within certain intervals. See custom_model documentation for
#' details.
#' 
#' @param hourtemps a list of two elements, with element 'hourtemps' being a
#' dataframe of hourly temperatures (e.g. produced by stack_hourly_temps). This
#' data frame must have a column for Year, a column for JDay (Julian date, or
#' day of the year), a column for Hour and a column for Temp (hourly
#' temperature). The second (optional) element is QC, which is a data.frame
#' indicating completeness of the dataset. This is automatically produced by
#' stack_hourly_temps.
#' @param Start_JDay the start date (in Julian date, or day of the year) of the
#' period, for which chill and heat should be quantified.
#' @param End_JDay the end date (in Julian date, or day of the year) of the
#' period, for which chill and heat should be quantified.
#' @param models named list of models that should be applied to the hourly
#' temperature data. These should be functions that take as input a vector of
#' hourly temperatures. This defaults to the set of models provided by the
#' chilling function.
#' @param misstolerance maximum percentage of values for a given season that
#' can be missing without the record being removed from the output. Defaults to
#' 50.
#' @return data frame showing totals for all specified models for the
#' respective periods for all seasons included in the temperature records.
#' Columns are Season, End_year (the year when the period ended) and Days (the
#' duration of the period), as well as one column per model, which receives the
#' same name as the function in the models list. If the weather input consisted
#' of a list with elements stack and QC, the output also contains columns from
#' QC that indicate the completeness of the weather record that the
#' calculations are based on.
#' @author Eike Luedeling
#' @references The chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords chill and heat calculation
#' @examples
#' 
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' 
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' df=data.frame(
#'   lower=c(-1000,1,2,3,4,5,6),
#'   upper=c(1,2,3,4,5,6,1000),
#'   weight=c(0,1,2,3,2,1,0))
#' 
#' custom<-function(x) step_model(x,df)
#' 
#' models<-list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,Chill_Portions=
#'   Dynamic_Model,GDH=GDH,custom=custom)
#' 
#' tempResponse(hourtemps,Start_JDay = 305,End_JDay = 60,models)
#' 
#' @export tempResponse
tempResponse <-
function (hourtemps,Start_JDay=1,End_JDay=366,models=list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,Chill_Portions=Dynamic_Model,GDH=GDH),misstolerance=50)             #hourtemps is a data frame with columns Year, JDay, Hour and Temp
     {
  if((length(names(hourtemps))==2) & ("hourtemps" %in% names(hourtemps))) {QC<-hourtemps$QC; hourtemps<-hourtemps$hourtemps} else QC<-NULL
  if(Start_JDay<End_JDay) {hourtemps[which(hourtemps$JDay>=Start_JDay&hourtemps$JDay<=End_JDay),"sea"]<-
                     hourtemps[which(hourtemps$JDay>=Start_JDay&hourtemps$JDay<=End_JDay),"Year"]} else
                        {hourtemps[which(hourtemps$JDay>=Start_JDay),"sea"]<-
                          hourtemps[which(hourtemps$JDay>=Start_JDay),"Year"]+1
                        hourtemps[which(hourtemps$JDay<=End_JDay),"sea"]<-
                          hourtemps[which(hourtemps$JDay<=End_JDay),"Year"]}
                         
      if(Start_JDay<End_JDay) {relevant_days<-Start_JDay:End_JDay} else
                              {relevant_days<-c(Start_JDay:366,1:End_JDay)}
      normal_lines<-which(!(hourtemps$JDay==Start_JDay&hourtemps$Hour==0))
      normal_lines<-normal_lines[which(normal_lines>1)]

      hourtemps<-hourtemps[which(!is.na(hourtemps[,"Temp"])),]

      
      for(m in 1:length(models))
        hourtemps[,names(models)[m]]<-do.call(models[[m]], list(hourtemps[,"Temp"]), quote = FALSE, envir = parent.frame())

      seasons<-unique(hourtemps$sea)
      seasons<-seasons[!is.na(seasons)]

#summarize all models
     chillout<-data.frame(Season=paste(seasons-1,"/",seasons,sep=""),End_year=seasons)
     
     if(End_JDay>=Start_JDay)
       dt<-sapply(seasons,function(x)
         difftime(ISOdate(x-1,12,31)+End_JDay*86400,ISOdate(x-1,12,31)+Start_JDay*86400))
     if(End_JDay<Start_JDay)
       dt<-sapply(seasons,function(x)
         difftime(ISOdate(x-1,12,31)+End_JDay*86400,ISOdate(x-2,12,31)+Start_JDay*86400))     
     chillout[,"Season_days"]<-dt+1
     chillout[,"Data_days"]<-sapply(seasons,function(x) length(which(hourtemps$sea==x))/24)

                                                              
for (sea in seasons)
    #if(sea==seasons[1])
     {#chillout<-data.frame(Season=paste(sea-1,"/",sea,sep=""),End_year=sea,Days=length(which(hourtemps$sea==sea))/24)
      seas<-hourtemps[which(hourtemps$sea==sea),]
      if("no_Tmin" %in% names(hourtemps)&"no_Tmax" %in% names(hourtemps))
        chillout[which(chillout$End_year==sea),"Interpolated_days"]<-length(which(seas$no_Tmin|seas$no_Tmax))/24
        chillout[,"Perc_complete"]<-NA
      if(sea==seasons[1]) last_end<-hourtemps[1,]-hourtemps[1,] else last_end<-hourtemps[min(which(hourtemps$sea==sea))-1,]
        for(m in 1:length(models))
          chillout[which(chillout$End_year==sea),names(models)[m]]<-seas[nrow(seas),names(models)[m]]-last_end[,names(models)[m]]
}
     if("no_Tmin" %in% names(hourtemps)&"no_Tmax" %in% names(hourtemps))
       chillout[,"Perc_complete"]<-(chillout[,"Data_days"]-chillout[,"Interpolated_days"])/chillout[,"Season_days"]*100 else
         chillout[,"Perc_complete"]<-chillout[,"Data_days"]/chillout[,"Season_days"]*100
           
     chillout<-chillout[which(chillout$Perc_complete>=100-misstolerance),]
     
     return(chillout)

}

