

Standard_summary=function (final_dat,bed_time,takeoff_time) {
  if(missing(takeoff_time)){
    final_dat=clean_time(final_dat)}else{
    final_dat=clean_time(final_dat,takeoff_time)
  }
if(missing(bed_time)){
table1=c()
# record.getup<-bed_time
# ###
# record.getup.time<-c()
# record.sleep.time<-c()
# for (kk in 1:nrow(record.getup))
# {
#   temp.getup.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,3] ),as.character(record.getup[kk,4]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
#   record.getup.time<-c(record.getup.time,temp.getup.time)
#   temp.sleep.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,5] ),as.character(record.getup[kk,6]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
#   record.sleep.time<-c(record.sleep.time,temp.sleep.time)
# }
  ###################################################
ll=0 #count for getup time
  for (i in unique(final_dat$new_date)){
    ll=ll+1
    final_dat_oneday=final_dat[final_dat$new_date==i,]
    # final_dat_oneday_clean<-final_dat_oneday[final_dat_oneday[,1]>record.getup.time[ll] &final_dat_oneday[,1]<record.sleep.time[ll], ]
final_dat_oneday_clean=final_dat_oneday
  time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))

  #hour_char<-as.numeric(format(as.POSIXlt(temp_mat$date_time*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0)),"%H"))
  ###
  temp_sed<-final_dat_oneday_clean[final_dat_oneday_clean$ActivityCode==0,]$Interval
  month<-as.numeric(format(time_char,"%m"))
   day<-as.numeric(format(time_char,"%d"))
   year<-as.numeric(format(time_char,"%Y"))
  hours_worn_total<- sum( final_dat_oneday_clean$Interval)/60/60
  #hours_awake<- (record_sleep_time[ll]-record_getup_time[ll])*24
  sed_hour<- sum(temp_sed) /60/60
  stand_hour<- sum(subset( final_dat_oneday_clean, final_dat_oneday_clean$ActivityCode==1)$Interval) /60/60
  step_hour<- sum(subset( final_dat_oneday_clean, final_dat_oneday_clean$ActivityCode==2)$Interval) /60/60
  ######## Code update here !
  step_count_total<- 2*nrow(subset( final_dat_oneday_clean, final_dat_oneday_clean$ActivityCode==2))  ##### in the event file, one row with ActivityCode==2 means two steps
  ########
 # gini_index<- gini(temp_sed)
  num_hour_over_3_METs<-  sum( final_dat_oneday_clean$Interval[( final_dat_oneday_clean$METs/ final_dat_oneday_clean$Interval)*60*60>3])/60/60
  MET_hours<- sum( final_dat_oneday_clean$METs)
  #valid_day<-ifelse( (hours_worn_total<10 & hours_worn_total/hours_awake>0_8 & step_count_total>200) | (hours_worn_total>=10 & step_count_total>200 ),1,0    )
  dayofweek<-as.numeric(format(time_char,"%w"))
  weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)
  #table1<-c(person,group_char,www,month,day,year,hours_worn_total,hours_awake,sed_hour,stand_hour,step_hour,num_changes_from_sed_to_non_sed,step_count_total,num_hour_over_3_METs,MET_hours,valid_day,dayofweek,weekday_or_weekend)
  perc_sedentary<- 100*sed_hour/hours_worn_total
  perc_stand<- 100*stand_hour/hours_worn_total
  perc_step<- 100*step_hour/hours_worn_total
  step_per_day<-step_count_total
  table<-cbind(year,month,day,dayofweek,weekday_or_weekend,hours_worn_total,sed_hour,stand_hour,step_hour,step_count_total,num_hour_over_3_METs,MET_hours,perc_sedentary,perc_stand,perc_step)

  table1=rbind(table1,table)

  }

 #colnames(table1)<-c("year","month","day","dayofweek","weekday_or_weekend","hours_worn_total","sed_hour","stand_hour","step_hour","step_count_total","num_hour_over_3_METs","MET_hours","perc_sedentary","perc_stand","perc_step")

 # table<-cbind(year,month,day,dayofweek,weekday_or_weekend,hours_worn_total,sed_hour,stand_hour,step_hour,num_changes_from_sed_to_non_sed,step_count_total,num_hour_over_3_METs,MET_hours)
  #colnames(table)<-c("year","month","day","dayofweek","weekday_or_weekend","hours_worn_total","sed_hour","stand_hour","step_hour","num_changes_from_sed_to_non_sed","step_count_total","num_hour_over_3_METs","MET_hours")

#    list(year=table1[,1],month=table1[,2], day=table1[,3],dayofweek=table1[,4], weekday_or_weekend=table1[,5],hours_worn_total=table1[,6], sed_hour=table1[,7], stand_hour=table1[,8], step_hour=table1[,9],step_count_total=table1[,10],num_hour_over_3_METs=table1[,11],MET_hours=table1[,12],perc_sedentary=table1[,13],perc_stand=table1[,14],perc_step=table1[,15])

  #table1_label<-c("id","group","week","month","day","year","hours_worn_total","hours_awake","sed_hour","stand_hour","step_hour","num_changes_from_sed_to_non_sed","step_count_total","num_hour_over_3_METs","MET_hours","valid_day","dayofweek","weekday_or_weekend")
  #num_hour_over_3_METs=num_hour_over_3_METs
  #return(table1)
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
    final_dat_oneday=final_dat[final_dat$new_date==i,]
    final_dat_oneday_clean<-final_dat_oneday[final_dat_oneday[,1]>record.getup.time[ll] &final_dat_oneday[,1]<record.sleep.time[ll], ]
        time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))

    #hour_char<-as.numeric(format(as.POSIXlt(temp_mat$date_time*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0)),"%H"))
    ###
    temp_sed<-final_dat_oneday_clean[final_dat_oneday_clean$ActivityCode==0,]$Interval
    month<-as.numeric(format(time_char,"%m"))
    day<-as.numeric(format(time_char,"%d"))
    year<-as.numeric(format(time_char,"%Y"))
    hours_worn_total<- sum( final_dat_oneday_clean$Interval)/60/60
    #hours_awake<- (record_sleep_time[ll]-record_getup_time[ll])*24
    sed_hour<- sum(temp_sed) /60/60
    stand_hour<- sum(subset( final_dat_oneday_clean, final_dat_oneday_clean$ActivityCode==1)$Interval) /60/60
    step_hour<- sum(subset( final_dat_oneday_clean, final_dat_oneday_clean$ActivityCode==2)$Interval) /60/60
    ######## Code update here !
    step_count_total<- 2*nrow(subset( final_dat_oneday_clean, final_dat_oneday_clean$ActivityCode==2))  ##### in the event file, one row with ActivityCode==2 means two steps
    ########
    # gini_index<- gini(temp_sed)
    num_hour_over_3_METs<-  sum( final_dat_oneday_clean$Interval[( final_dat_oneday_clean$METs/ final_dat_oneday_clean$Interval)*60*60>3])/60/60
    MET_hours<- sum( final_dat_oneday_clean$METs)
    #valid_day<-ifelse( (hours_worn_total<10 & hours_worn_total/hours_awake>0_8 & step_count_total>200) | (hours_worn_total>=10 & step_count_total>200 ),1,0    )
    dayofweek<-as.numeric(format(time_char,"%w"))
    weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)
    #table1<-c(person,group_char,www,month,day,year,hours_worn_total,hours_awake,sed_hour,stand_hour,step_hour,num_changes_from_sed_to_non_sed,step_count_total,num_hour_over_3_METs,MET_hours,valid_day,dayofweek,weekday_or_weekend)
    perc_sedentary<- 100*sed_hour/hours_worn_total
    perc_stand<- 100*stand_hour/hours_worn_total
    perc_step<- 100*step_hour/hours_worn_total
    step_per_day<-step_count_total
    table<-cbind(year,month,day,dayofweek,weekday_or_weekend,hours_worn_total,sed_hour,stand_hour,step_hour,step_count_total,num_hour_over_3_METs,MET_hours,perc_sedentary,perc_stand,perc_step)

    table1=rbind(table1,table)

  }

  colnames(table1)<-c("year","month","day","dayofweek","weekday_or_weekend","hours_worn_total","sed_hour","stand_hour","step_hour","step_count_total","num_hour_over_3_METs","MET_hours","perc_sedentary","perc_stand","perc_step")

  # table<-cbind(year,month,day,dayofweek,weekday_or_weekend,hours_worn_total,sed_hour,stand_hour,step_hour,num_changes_from_sed_to_non_sed,step_count_total,num_hour_over_3_METs,MET_hours)
  #colnames(table)<-c("year","month","day","dayofweek","weekday_or_weekend","hours_worn_total","sed_hour","stand_hour","step_hour","num_changes_from_sed_to_non_sed","step_count_total","num_hour_over_3_METs","MET_hours")
}
   list(year=table1[,1],month=table1[,2], day=table1[,3],dayofweek=table1[,4], weekday_or_weekend=table1[,5],hours_worn_total=table1[,6], sed_hour=table1[,7], stand_hour=table1[,8], step_hour=table1[,9],step_count_total=table1[,10],num_hour_over_3_METs=table1[,11],MET_hours=table1[,12],perc_sedentary=table1[,13],perc_stand=table1[,14],perc_step=table1[,15])

  #table1_label<-c("id","group","week","month","day","year","hours_worn_total","hours_awake","sed_hour","stand_hour","step_hour","num_changes_from_sed_to_non_sed","step_count_total","num_hour_over_3_METs","MET_hours","valid_day","dayofweek","weekday_or_weekend")
  #num_hour_over_3_METs=num_hour_over_3_METs
  #return(table1)

}
#' Standard Physical Activity Summary
#'
#' Summarize standard activity measures
#' @param final_dat Raw event file, will be cleaned in this function. Event file is required for this function.
#' @param takeoff_time Take on and off time log, reported by participants. Log is not required for this function.
#' @param bed_time Sleep and wake up time log, reported by participants. Log is not required for this function.
#' @return \code{Year} The calendar year of recorded event
#' @return \code{Month} The calendar month of recorded event
#' @return \code{Day}   The calendar day of recorded event
#' @return \code{Dayofweek} The day of that week
#' @return \code{Weekday_or_weekend}  The recored event date is a weekday or weekend (0 for weekday, 1 for weekend)
#' @return \code{hours_worn_total} Total worn hours
#' @return \code{sed_hour}  Total sedantary hours
#' @return \code{stand_hour} Total stand hours
#' @return \code{step_hour} Total step hours
#' @return \code{step_count_total} Total Step count
#' @return \code{num_hour_over_3_METs} Number of hours that Metabolic Equivalent of Task (MET) is over 3
#' @return \code{MET_hours}  Total METs
#' @return \code{perc_sedentary} Proportion of sedentary
#' @return \code{perc_stand} Proportion of standing
#' @return \code{perc_step} Proportion of stepping
#' @details All numbers are calculated in the given time period (day, hour, etc.). Total sedentary/standing/stepping hours are obtained from the summation of the duration times for sedentary/standing/stepping activities in the given time period.
#' @examples
#' #For CRAN less than 5s running time policy, we only select the first day to run.
#' r1=Standard(sample_event[1:3095,],sample_bed_time[1,],sample_takeon_log[1,])
#' summary(r1)
#' @export
#'
Standard=function(final_dat,bed_time,takeoff_time) UseMethod("Standard")
#' @export

Standard.default=function(final_dat,bed_time,takeoff_time)
{
  out=Standard_summary(final_dat,bed_time,takeoff_time)
  out$call=match.call()
  class(out)="Standard"
  out
}
#' @export

print.Standard=function(x,...)
{cat("Call:\n")
  print(x$call)
}
#' @export


summary.Standard=function(object,...)
{
  TAB=cbind(object$year,object$month,object$day,object$dayofweek,object$weekday_or_weekend,object$hours_worn_total,object$sed_hour,object$stand_hour,object$step_hour,object$step_count_total,object$num_hour_over_3_METs,object$MET_hours,object$perc_sedentary,object$perc_stand,object$perc_step)

  colnames(TAB)=c("Year","Month","Day","Dayofweek","Weekday_or_weekend","hours_worn_total","sed_hour","stand_hour","step_hour","step_count_total","num_hour_over_3_METs","MET_hours","perc_sedentary","perc_stand","perc_step")
  res <- list(call=object$call,
              Table=TAB)
  class(res) <- "summary.Standard"
  res
}

