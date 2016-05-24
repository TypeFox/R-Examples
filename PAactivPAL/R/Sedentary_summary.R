
Sedentary_summary=function (final_dat,bed_time,takeoff_time) {
  options(digits = 10)

  if(missing(takeoff_time)){
    final_dat=clean_time(final_dat)}else{
      final_dat=clean_time(final_dat,takeoff_time)
    }
  if(missing(bed_time)){
    table1=c()
#     record.getup<-bed_time
#     ###
#     record.getup.time<-c()
#     record.sleep.time<-c()
#     for (kk in 1:nrow(record.getup))
#     {
#       temp.getup.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,3] ),as.character(record.getup[kk,4]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
#       record.getup.time<-c(record.getup.time,temp.getup.time)
#       temp.sleep.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,5] ),as.character(record.getup[kk,6]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
#       record.sleep.time<-c(record.sleep.time,temp.sleep.time)
#     }
    ###################################################

    for (i in unique(final_dat$new_date)){

      temp_mat_oneday=final_dat[final_dat$new_date==i,]
      # temp_mat_oneday<-temp_mat_oneday[temp_mat_oneday[,1]>record.getup.time[ll] &temp_mat_oneday[,1]<record.sleep.time[ll], ]

      time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))
      month<-as.numeric(format(time_char,"%m"))
      day<-as.numeric(format(time_char,"%d"))
      year<-as.numeric(format(time_char,"%Y"))
      dayofweek<-as.numeric(format(time_char,"%w"))
      weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)
      hour_char<-as.numeric(format(as.POSIXlt(temp_mat_oneday$date_time*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0)),"%H"))
      temp_sed<-subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0)$Interval
      length_temp_sed<-length(temp_sed)
      num_changes_from_sed_to_non_sed<- length_temp_sed

      sed_hour<- sum(temp_sed) /60/60
      length_temp_sed<-length(temp_sed)
      num_changes_from_sed_to_non_sed<- length_temp_sed
      total_sed_time<-sed_hour
      #total_number_of_sed_bouts<-num_changes_from_sed_to_non_sed
      mean_sed_bout_length<- mean(temp_sed) /60/60
      prop_of_sed_time_greater_20min<- 100*length(temp_sed[temp_sed>20*60])/length_temp_sed
      prop_of_sed_time_greater_60min<- 100*length(temp_sed[temp_sed>60*60])/length_temp_sed
      prop_of_sed_time_greater_120min<- 100*length(temp_sed[temp_sed>120*60])/length_temp_sed

      total_sed_time_greater_20min<- sum(temp_sed[temp_sed>20*60])/60/60
      total_sed_time_greater_60min<- sum(temp_sed[temp_sed>60*60])/60/60
      total_sed_time_greater_120min<- sum(temp_sed[temp_sed>120*60])/60/60

      sed_time_in_30_and_60=sum(temp_sed[temp_sed>=30*60 & temp_sed<=60*60])/60
      sed_num_in_30_and_60=length(temp_sed[temp_sed>=30*60 & temp_sed<=60*60])
      quantile_temp<-quantile(temp_sed, probs = c(0.05,0.25,0.5,0.75,0.95))/60/60
      percentile_sed_time_5<- quantile_temp[1]
      percentile_sed_time_25<- quantile_temp[2]
      percentile_sed_time_50<- quantile_temp[3]
      percentile_sed_time_75<- quantile_temp[4]
      percentile_sed_time_95<- quantile_temp[5]

      alpha_sed<- 1+ 1/mean(log(temp_sed/ min(temp_sed)))
      #gini_index_sed<- gini(temp_sed)

      prop_sed_time_6_12<- 100*sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0 & hour_char>=6 & hour_char<12)$Interval) /(sum(subset(temp_mat_oneday,hour_char>=6 & hour_char<12)$Interval)+0.0001)   ###prevent this value is zero
      prop_sed_time_12_18<- 100*sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0 & hour_char>=12 & hour_char<18)$Interval) /(sum(subset(temp_mat_oneday,hour_char>=12 & hour_char<18)$Interval)+0.0001)
      prop_sed_time_18_22<- 100*sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0 & hour_char>=18 & hour_char<22)$Interval) /(sum(subset(temp_mat_oneday,hour_char>=18 & hour_char<22)$Interval)+0.0001)

      break_per_day<-num_changes_from_sed_to_non_sed
      break_rate<-break_per_day/sed_hour
      table<- cbind(year,month,day,dayofweek,weekday_or_weekend,break_per_day,break_rate,mean_sed_bout_length,prop_of_sed_time_greater_20min,prop_of_sed_time_greater_60min,prop_of_sed_time_greater_120min,total_sed_time_greater_20min,total_sed_time_greater_60min,total_sed_time_greater_120min,percentile_sed_time_5,percentile_sed_time_25,percentile_sed_time_50,percentile_sed_time_75,percentile_sed_time_95,alpha_sed,prop_sed_time_6_12,prop_sed_time_12_18,prop_sed_time_18_22,sed_time_in_30_and_60,sed_num_in_30_and_60)
      row.names(table)=NULL
      table1=rbind(table1,table)
      #print(table1)
    }
  }else{
    table1=c()
  record.getup<-bed_time
  ###
  record.getup.time<-c()
  record.sleep.time<-c()
  for (kk in 1:nrow(record.getup))
  {
    temp.getup.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,3] ),as.character(record.getup[kk,4]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
    record.getup.time<-c(record.getup.time,temp.getup.time)
    temp.sleep.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,5] ),as.character(record.getup[kk,6]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
    record.sleep.time<-c(record.sleep.time,temp.sleep.time)
  }
  ###################################################
  ll=0 #count for getup time

  for (i in unique(final_dat$new_date)){
    ll=ll+1
    temp_mat_oneday=final_dat[final_dat$new_date==i,]
    temp_mat_oneday<-temp_mat_oneday[temp_mat_oneday[,1]>record.getup.time[ll] &temp_mat_oneday[,1]<record.sleep.time[ll], ]

    time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))
    month<-as.numeric(format(time_char,"%m"))
    day<-as.numeric(format(time_char,"%d"))
    year<-as.numeric(format(time_char,"%Y"))
    dayofweek<-as.numeric(format(time_char,"%w"))
    weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)
  hour_char<-as.numeric(format(as.POSIXlt(temp_mat_oneday$date_time*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0)),"%H"))
  temp_sed<-subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0)$Interval
  length_temp_sed<-length(temp_sed)
  num_changes_from_sed_to_non_sed<- length_temp_sed

  sed_hour<- sum(temp_sed) /60/60
  length_temp_sed<-length(temp_sed)
  num_changes_from_sed_to_non_sed<- length_temp_sed
  total_sed_time<-sed_hour
  #total_number_of_sed_bouts<-num_changes_from_sed_to_non_sed
  mean_sed_bout_length<- mean(temp_sed) /60/60
  prop_of_sed_time_greater_20min<- 100*length(temp_sed[temp_sed>20*60])/length_temp_sed
  prop_of_sed_time_greater_60min<- 100*length(temp_sed[temp_sed>60*60])/length_temp_sed
  prop_of_sed_time_greater_120min<- 100*length(temp_sed[temp_sed>120*60])/length_temp_sed

  total_sed_time_greater_20min<- sum(temp_sed[temp_sed>20*60])/60/60
  total_sed_time_greater_60min<- sum(temp_sed[temp_sed>60*60])/60/60
  total_sed_time_greater_120min<- sum(temp_sed[temp_sed>120*60])/60/60

  sed_time_in_30_and_60=sum(temp_sed[temp_sed>=30*60 & temp_sed<=60*60])/60
  sed_num_in_30_and_60=length(temp_sed[temp_sed>=30*60 & temp_sed<=60*60])
  quantile_temp<-quantile(temp_sed, probs = c(0.05,0.25,0.5,0.75,0.95))/60/60
  percentile_sed_time_5<- quantile_temp[1]
  percentile_sed_time_25<- quantile_temp[2]
  percentile_sed_time_50<- quantile_temp[3]
  percentile_sed_time_75<- quantile_temp[4]
  percentile_sed_time_95<- quantile_temp[5]

  alpha_sed<- 1+ 1/mean(log(temp_sed/ min(temp_sed)))
  #gini_index_sed<- gini(temp_sed)

  prop_sed_time_6_12<- 100*sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0 & hour_char>=6 & hour_char<12)$Interval) /(sum(subset(temp_mat_oneday,hour_char>=6 & hour_char<12)$Interval)+0.0001)   ###prevent this value is zero
  prop_sed_time_12_18<- 100*sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0 & hour_char>=12 & hour_char<18)$Interval) /(sum(subset(temp_mat_oneday,hour_char>=12 & hour_char<18)$Interval)+0.0001)
  prop_sed_time_18_22<- 100*sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==0 & hour_char>=18 & hour_char<22)$Interval) /(sum(subset(temp_mat_oneday,hour_char>=18 & hour_char<22)$Interval)+0.0001)

  break_per_day<-num_changes_from_sed_to_non_sed
  break_rate<-break_per_day/sed_hour
  table<- cbind(year,month,day,dayofweek,weekday_or_weekend,break_per_day,break_rate,mean_sed_bout_length,prop_of_sed_time_greater_20min,prop_of_sed_time_greater_60min,prop_of_sed_time_greater_120min,total_sed_time_greater_20min,total_sed_time_greater_60min,total_sed_time_greater_120min,percentile_sed_time_5,percentile_sed_time_25,percentile_sed_time_50,percentile_sed_time_75,percentile_sed_time_95,alpha_sed,prop_sed_time_6_12,prop_sed_time_12_18,prop_sed_time_18_22,sed_time_in_30_and_60,sed_num_in_30_and_60)
  row.names(table)=NULL
  table1=rbind(table1,table)
  #print(table1)

  }
  }

  #colnames(table1)<- c("year","month","day","dayofweek","weekday_or_weekend","break_per_day","break_rate","mean_sed_bout_length" ,"prop_of_sed_time_greater_20min","prop_of_sed_time_greater_60min","prop_of_sed_time_greater_120min","total_sed_time_greater_20min","total_sed_time_greater_60min","total_sed_time_greater_120min","percentile_sed_time_5","percentile_sed_time_25","percentile_sed_time_50","percentile_sed_time_75","percentile_sed_time_95","alpha_sed","prop_sed_time_6_12","prop_sed_time_12_18","prop_sed_time_18_22","sed_time_in_30_and_60","sed_num_in_30_and_60")
  list(year=table1[,1],month=table1[,2], day=table1[,3],dayofweek=table1[,4],weekday_or_weekend=table1[,5], break_per_day=table1[,6],break_rate=table1[,7], mean_sed_bout_length=table1[,8], prop_of_sed_time_greater_20min=table1[,9], prop_of_sed_time_greater_60min=table1[,10],prop_of_sed_time_greater_120min=table1[,11],total_sed_time_greater_20min=table1[,12],total_sed_time_greater_60min=table1[,13],total_sed_time_greater_120min=table1[,14],percentile_sed_time_5=table1[,15],percentile_sed_time_25=table1[,16],percentile_sed_time_50=table1[,17],percentile_sed_time_75=table1[,18],percentile_sed_time_95=table1[,19],alpha_sed=table1[,20],prop_sed_time_6_12=table1[,21],prop_sed_time_12_18=table1[,22],prop_sed_time_18_22=table1[,23],sed_time_in_30_and_60=table1[,24],sed_num_in_30_and_60=table1[,25])

#   out=list( total_number_of_activity_bouts=total_number_of_activity_bouts,mean_activity_bout_length=mean_activity_bout_length,prop_of_activity_time_greater_5min=prop_of_activity_time_greater_5min,prop_of_activity_time_greater_10min=prop_of_activity_time_greater_10min,prop_of_activity_time_greater_30min=prop_of_activity_time_greater_30min,total_activity_time_greater_5min=total_activity_time_greater_5min,total_activity_time_greater_10min=total_activity_time_greater_10min,total_activity_time_greater_30min=total_activity_time_greater_30min,percentile_activity_time_5=percentile_activity_time_5,percentile_activity_time_25=percentile_activity_time_25,percentile_activity_time_50=percentile_activity_time_50,percentile_activity_time_75=percentile_activity_time_75,percentile_activity_time_95=percentile_activity_time_95,alpha_activity=alpha_activity,stepping_to_standing_ratio=stepping_to_standing_ratio,table=table)
  #return(table1)
}

#' Sedentary Physical Activity Summary
#'
#' Summarize sedentary features using proportions and percentiles
#' @param final_dat Raw event file, will be cleaned in this function. Event file is required for this function.
#' @param takeoff_time Take on and off time log, reported by participants. Log is not required for this function.
#' @param bed_time Sleep and wake up time log, reported by participants. Log is not required for this function.
#' @return \code{Year} The calendar year of recorded event
#' @return \code{Month} The calendar month of recorded event
#' @return \code{Day}   The calendar day of recorded event
#' @return \code{Dayofweek} The day of that week
#' @return \code{Weekday_or_weekend}  The recored event date is a weekday or weekend (0 for weekday, 1 for weekend)
#' @return \code{break_per_day} Number of interrupting sedentary behavior during the given period of time
#' @return \code{break_rate} Rate of interrupting sedentary behavior
#' @return \code{mean_sed_bout_length} Mean sedentary bout length. It is the average of the duration time for all sedentary activities
#' @return \code{prop_of_sed_time_greater_20min} Proportions of sedentary bout greater than 20 minutes
#' @return \code{prop_of_sed_time_greater_60min} Proportions of sedentary bout greater than 60 minutes
#' @return \code{prop_of_sed_time_greater_120min} Proportions of sedentary bout greater than 120 minutes
#' @return \code{total_sed_time_greater_20min} Total sedentary time greater than 20 minutes
#' @return \code{total_sed_time_greater_60min} Total sedentary time greater than 60 minutes
#' @return \code{total_sed_time_greater_120min} Total sedentary time greater than 120 minutes
#' @return \code{percentile_sed_time_5}  5\% Percentile of sedentary time
#' @return \code{percentile_sed_time_25} 25\% Percentile of sedentary time
#' @return \code{percentile_sed_time_50} 50\% Percentile of sedentary time
#' @return \code{percentile_sed_time_75} 75\% Percentile of sedentary time
#' @return \code{percentile_sed_time_95} 95\% Percentile of sedentary time
#' @return \code{alpha_sed} alpha of sedentary time, see details
#' @return \code{prop_sed_time_6_12} Proportions of sedentary time between 6:00-12:00
#' @return \code{prop_sed_time_12_18} Proportions of sedentary time between 12:00-18:00
#' @return \code{prop_sed_time_18_22} Proportions of sedentary time between 18:00-22:00
#' @return \code{sed_time_in_30_and_60} Minutes Spent in Sitting/Lying Bouts (at least 30 and less 60 minutes in duration)
#' @return \code{sed_num_in_30_and_60} Number of in Sitting/Lying Bouts (at least 30 and less 60 minutes in duration)
#' @details Proportions of sedentary bout greater than 20/60/120 minutes is the ratio of the number of sedentary bouts greater than 20/60/120 minutes to the total number of sedentary recordings.
#' @details Total sedentary time greater than 20/60/120 minutes is the summation of the sedentary durations which are greater than 20/60/120 minutes.
#' @details To calculate 5\%/25\%/75\%/95\% percentile of sedentary time, all of the recorded sedentary durations are listed and R function \emph{quantile} is used to find the percentiles.
#' @details  alpha_sed is defined by \code{1+1/M}, where \code{M} is the average of \code{log(sedentary bout length /minimum sedentary bout length)}.
#' @details Proportions of sedentary time between 6:00-12:00/12:00-18:00/18:00-22:00 is the ratio of the sedentary durations to the total activity durations between 6:00-12:00/12:00-18:00/18:00-22:00.
#' @importFrom stats quantile
#' @examples
#' #For CRAN less than 5s running time policy, we only select the first day to run.
#' r2=Sedentary(sample_event[1:3095,],sample_bed_time[1,],sample_takeon_log[1,])
#' summary(r2)
#' @export
#'

Sedentary=function(final_dat,bed_time,takeoff_time) UseMethod("Sedentary")
#' @export

Sedentary.default=function(final_dat,bed_time,takeoff_time)
{
  out=Sedentary_summary(final_dat,bed_time,takeoff_time)
  out$call=match.call()
  class(out)="Sedentary"
  out
}
#' @export

print.Sedentary=function(x,...)
{cat("Call:\n")
  print(x$call)
}

#' @export

summary.Sedentary=function(object,...)
{
  TAB=cbind(object$year,object$month,object$day,object$dayofweek,object$weekday_or_weekend,object$break_per_day,object$break_rate,object$mean_sed_bout_length,object$prop_of_sed_time_greater_20min,object$prop_of_sed_time_greater_60min,object$prop_of_sed_time_greater_120min,object$total_sed_time_greater_20min,object$total_sed_time_greater_60min,object$total_sed_time_greater_120min,object$percentile_sed_time_5,object$percentile_sed_time_25,object$percentile_sed_time_50,object$percentile_sed_time_75,object$percentile_sed_time_95,object$alpha_sed,object$prop_sed_time_6_12,object$prop_sed_time_12_18,object$prop_sed_time_18_22,object$sed_time_in_30_and_60,object$sed_num_in_30_and_60)

  colnames(TAB)=c("Year","Month","Day","Dayofweek","Weekday_or_weekend","break_per_day","break_rate","mean_sed_bout_length" ,"prop_of_sed_time_greater_20min","prop_of_sed_time_greater_60min","prop_of_sed_time_greater_120min","total_sed_time_greater_20min","total_sed_time_greater_60min","total_sed_time_greater_120min","percentile_sed_time_5","percentile_sed_time_25","percentile_sed_time_50","percentile_sed_time_75","percentile_sed_time_95","alpha_sed","prop_sed_time_6_12","prop_sed_time_12_18","prop_sed_time_18_22","sed_time_in_30_and_60","sed_num_in_30_and_60")
  res <- list(call=object$call,
              Table=TAB)
  class(res) <- "summary.Sedentary"
  res
}



