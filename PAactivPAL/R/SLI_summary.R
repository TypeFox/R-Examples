

SLI_summary=function (final_dat,bed_time,takeoff_time) {
  if(missing(takeoff_time)){
    final_dat=clean_time(final_dat)}else{
      final_dat=clean_time(final_dat,takeoff_time)
    }
  if(missing(bed_time)){
    table1=c()

#   record.getup<-bed_time
#   ###
#   record.getup.time<-c()
#   record.sleep.time<-c()
#   for (kk in 1:nrow(record.getup))
#   {
#     temp.getup.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,3] ),as.character(record.getup[kk,4]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
#     record.getup.time<-c(record.getup.time,temp.getup.time)
#     temp.sleep.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record.getup[kk,5] ),as.character(record.getup[kk,6]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
#     record.sleep.time<-c(record.sleep.time,temp.sleep.time)
#   }

  temp_mat_for_activity<- final_dat
  temp_mat_for_activity$Activity[temp_mat_for_activity$Activity==1]<-2
  end_pos<-cumsum(rle(temp_mat_for_activity$Activity)$lengths)
  start_pos<-c(0,end_pos[1:(length(end_pos)-1)])+1
  ############### for each runs, handle the data
  handle_runs<- sapply(1:length(end_pos),function(x,data_mat=temp_mat_for_activity)
  {
    select_data<-data_mat[start_pos[x]:end_pos[x],]
    combine_data<- c(min(select_data$date_time),sum(select_data$Interval),max(select_data$Activity),sum(select_data$METs),min(select_data$new_date) )
    return(combine_data)
  }, simplify=F
  )
  ############### combine each run
  combined_temp_mat_for_activity<-data.frame(do.call(rbind,handle_runs))
  colnames(combined_temp_mat_for_activity)<-c("date_time", "Interval", "ActivityCode","METs","new_date")
  ###############
  table1=c()
  ###################################################
  ll=0 #count for getup time

  for (i in unique(combined_temp_mat_for_activity$new_date)){
    ll=ll+1
    temp_mat_oneday=final_dat[final_dat$new_date==i,]
#     temp_mat_oneday<-temp_mat_oneday[temp_mat_oneday[,1]>record.getup.time[ll] &temp_mat_oneday[,1]<record.sleep.time[ll], ]

    final_dat_oneday=combined_temp_mat_for_activity[combined_temp_mat_for_activity$new_date==i,]
   # final_dat_oneday<-final_dat_oneday[final_dat_oneday[,1]>record.getup.time[ll] &final_dat_oneday[,1]<record.sleep.time[ll], ]

    time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))
    month<-as.numeric(format(time_char,"%m"))
    day<-as.numeric(format(time_char,"%d"))
    year<-as.numeric(format(time_char,"%Y"))
    dayofweek<-as.numeric(format(time_char,"%w"))
    weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)


    ###############Calculation
    ###############
    ###
    temp_activity<-subset(final_dat_oneday,final_dat_oneday$ActivityCode==2)$Interval
    ###
    length_temp_activity<-length(temp_activity)
    total_number_of_activity_bouts<- length_temp_activity
    mean_activity_bout_length<- mean(temp_activity) /60/60

    prop_of_activity_time_greater_5min<- 100*length(temp_activity[temp_activity>5*60])/length_temp_activity
    prop_of_activity_time_greater_10min<- 100*length(temp_activity[temp_activity>10*60])/length_temp_activity
    prop_of_activity_time_greater_30min<- 100*length(temp_activity[temp_activity>30*60])/length_temp_activity

    total_activity_time_greater_5min<- sum(temp_activity[temp_activity>5*60])/60/60
    total_activity_time_greater_10min<- sum(temp_activity[temp_activity>10*60])/60/60
    total_activity_time_greater_30min<- sum(temp_activity[temp_activity>30*60])/60/60

    quantile_activity_temp<-quantile(temp_activity, probs = c(0.05,0.25,0.5,0.75,0.95))/60/60
    percentile_activity_time_5<- quantile_activity_temp[1]
    percentile_activity_time_25<- quantile_activity_temp[2]
    percentile_activity_time_50<- quantile_activity_temp[3]
    percentile_activity_time_75<- quantile_activity_temp[4]
    percentile_activity_time_95<- quantile_activity_temp[5]

    alpha_activity<- 1+ 1/mean(log(temp_activity/ min(temp_activity)))
  #  gini_index_activity<- gini(temp_activity)
    step_hour<- sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==2)$Interval) /60/60
    stand_hour<- sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==1)$Interval) /60/60

    stepping_to_standing_ratio<- step_hour/stand_hour

    table<- cbind(year,month,day,dayofweek,weekday_or_weekend,total_number_of_activity_bouts,mean_activity_bout_length,prop_of_activity_time_greater_5min,prop_of_activity_time_greater_10min,prop_of_activity_time_greater_30min,total_activity_time_greater_5min,total_activity_time_greater_10min,total_activity_time_greater_30min,percentile_activity_time_5,percentile_activity_time_25,percentile_activity_time_50,percentile_activity_time_75,percentile_activity_time_95,alpha_activity,stepping_to_standing_ratio)
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

    temp_mat_for_activity<- final_dat
    temp_mat_for_activity$Activity[temp_mat_for_activity$Activity==1]<-2
    end_pos<-cumsum(rle(temp_mat_for_activity$Activity)$lengths)
    start_pos<-c(0,end_pos[1:(length(end_pos)-1)])+1
    ############### for each runs, handle the data
    handle_runs<- sapply(1:length(end_pos),function(x,data_mat=temp_mat_for_activity)
    {
      select_data<-data_mat[start_pos[x]:end_pos[x],]
      combine_data<- c(min(select_data$date_time),sum(select_data$Interval),max(select_data$Activity),sum(select_data$METs),min(select_data$new_date) )
      return(combine_data)
    }, simplify=F
    )
    ############### combine each run
    combined_temp_mat_for_activity<-data.frame(do.call(rbind,handle_runs))
    colnames(combined_temp_mat_for_activity)<-c("date_time", "Interval", "ActivityCode","METs","new_date")
    ###############
    table1=c()
    ###################################################
    ll=0 #count for getup time

    for (i in unique(combined_temp_mat_for_activity$new_date)){
      ll=ll+1
      temp_mat_oneday=final_dat[final_dat$new_date==i,]
          temp_mat_oneday<-temp_mat_oneday[temp_mat_oneday[,1]>record.getup.time[ll] &temp_mat_oneday[,1]<record.sleep.time[ll], ]

      final_dat_oneday=combined_temp_mat_for_activity[combined_temp_mat_for_activity$new_date==i,]
      final_dat_oneday<-final_dat_oneday[final_dat_oneday[,1]>record.getup.time[ll] &final_dat_oneday[,1]<record.sleep.time[ll], ]

      time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))
      month<-as.numeric(format(time_char,"%m"))
      day<-as.numeric(format(time_char,"%d"))
      year<-as.numeric(format(time_char,"%Y"))
      dayofweek<-as.numeric(format(time_char,"%w"))
      weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)


      ###############Calculation
      ###############
      ###
      temp_activity<-subset(final_dat_oneday,final_dat_oneday$ActivityCode==2)$Interval
      ###
      length_temp_activity<-length(temp_activity)
      total_number_of_activity_bouts<- length_temp_activity
      mean_activity_bout_length<- mean(temp_activity) /60/60

      prop_of_activity_time_greater_5min<- 100*length(temp_activity[temp_activity>5*60])/length_temp_activity
      prop_of_activity_time_greater_10min<- 100*length(temp_activity[temp_activity>10*60])/length_temp_activity
      prop_of_activity_time_greater_30min<- 100*length(temp_activity[temp_activity>30*60])/length_temp_activity

      total_activity_time_greater_5min<- sum(temp_activity[temp_activity>5*60])/60/60
      total_activity_time_greater_10min<- sum(temp_activity[temp_activity>10*60])/60/60
      total_activity_time_greater_30min<- sum(temp_activity[temp_activity>30*60])/60/60

      quantile_activity_temp<-quantile(temp_activity, probs = c(0.05,0.25,0.5,0.75,0.95))/60/60
      percentile_activity_time_5<- quantile_activity_temp[1]
      percentile_activity_time_25<- quantile_activity_temp[2]
      percentile_activity_time_50<- quantile_activity_temp[3]
      percentile_activity_time_75<- quantile_activity_temp[4]
      percentile_activity_time_95<- quantile_activity_temp[5]

      alpha_activity<- 1+ 1/mean(log(temp_activity/ min(temp_activity)))
      #  gini_index_activity<- gini(temp_activity)
      step_hour<- sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==2)$Interval) /60/60
      stand_hour<- sum(subset(temp_mat_oneday,temp_mat_oneday$ActivityCode==1)$Interval) /60/60

      stepping_to_standing_ratio<- step_hour/stand_hour

      table<- cbind(year,month,day,dayofweek,weekday_or_weekend,total_number_of_activity_bouts,mean_activity_bout_length,prop_of_activity_time_greater_5min,prop_of_activity_time_greater_10min,prop_of_activity_time_greater_30min,total_activity_time_greater_5min,total_activity_time_greater_10min,total_activity_time_greater_30min,percentile_activity_time_5,percentile_activity_time_25,percentile_activity_time_50,percentile_activity_time_75,percentile_activity_time_95,alpha_activity,stepping_to_standing_ratio)
      row.names(table)=NULL
      table1=rbind(table1,table)
      #print(table1)

    }

#   colnames(table1)<- c("year","month","day","dayofweek","weekday_or_weekend","total_number_of_activity_bouts","mean_activity_bout_length","prop_of_activity_time_greater_5min","prop_of_activity_time_greater_10min","prop_of_activity_time_greater_30min","total_activity_time_greater_5min","total_activity_time_greater_10min","total_activity_time_greater_30min","percentile_activity_time_5","percentile_activity_time_25","percentile_activity_time_50","percentile_activity_time_75","percentile_activity_time_95","alpha_activity","stepping_to_standing_ratio")

  #out=list( total_number_of_activity_bouts=total_number_of_activity_bouts,mean_activity_bout_length=mean_activity_bout_length,prop_of_activity_time_greater_5min=prop_of_activity_time_greater_5min,prop_of_activity_time_greater_10min=prop_of_activity_time_greater_10min,prop_of_activity_time_greater_30min=prop_of_activity_time_greater_30min,total_activity_time_greater_5min=total_activity_time_greater_5min,total_activity_time_greater_10min=total_activity_time_greater_10min,total_activity_time_greater_30min=total_activity_time_greater_30min,percentile_activity_time_5=percentile_activity_time_5,percentile_activity_time_25=percentile_activity_time_25,percentile_activity_time_50=percentile_activity_time_50,percentile_activity_time_75=percentile_activity_time_75,percentile_activity_time_95=percentile_activity_time_95,alpha_activity=alpha_activity,stepping_to_standing_ratio=stepping_to_standing_ratio,table=table)
  }
  list(year=table1[,1],month=table1[,2], day=table1[,3],dayofweek=table1[,4], weekday_or_weekend=table1[,5],total_number_of_activity_bouts=table1[,6], mean_activity_bout_length=table1[,7],prop_of_activity_time_greater_5min=table1[,8], prop_of_activity_time_greater_10min=table1[,9],prop_of_activity_time_greater_30min=table1[,10],total_activity_time_greater_5min=table1[,11],total_activity_time_greater_10min=table1[,12],total_activity_time_greater_30min=table1[,13],percentile_activity_time_5=table1[,14],percentile_activity_time_25=table1[,15],percentile_activity_time_50=table1[,16],percentile_activity_time_75=table1[,17],percentile_activity_time_95=table1[,18],alpha_activity=table1[,19],stepping_to_standing_ratio=table1[,20])

}
#' Standing/Light Intensity Physical Activity Summary
#'
#' Summarize standing and light intensity activity measures
#' @param final_dat Raw event file, will be cleaned in this function. Event file is required for this function.
#' @param takeoff_time Take on and off time log, reported by participants. Log is not required for this function.
#' @param bed_time Sleep and wake up time log, reported by participants. Log is not required for this function.
#' @return \code{Year} The calendar year of recorded event
#' @return \code{Month} The calendar month of recorded event
#' @return \code{Day}   The calendar day of recorded event
#' @return \code{Dayofweek} The day of that week
#' @return \code{Weekday_or_weekend}  The recored event date is a weekday or weekend (0 for weekday, 1 for weekend)
#' @return \code{total_number_of_activity_bouts} Total events that are standing or stepping
#' @return \code{total_number_of_activity_bouts} Total number of events that are standing or stepping
#' @return \code{mean_activity_bout_length} Average length of activity bout
#' @return \code{prop_of_activity_time_greater_5min} The ratio of the number of active bouts greater than 5 min to the total number of active events
#' @return \code{prop_of_activity_time_greater_10min} The ratio of the number of active bouts greater than 10 min to the total number of active events.
#' @return \code{prop_of_activity_time_greater_30min} The ratio of the number of active bouts greater than 30 min to the total number of active events.
#' @return \code{total_activity_time_greater_5min} The summation of the active durations which are greater than 5 min
#' @return \code{total_activity_time_greater_10min} The summation of the active durations which are greater than 10 min
#' @return \code{total_activity_time_greater_30min} The summation of the active durations which are greater than 30 min
#' @return \code{percentile_activity_time_5} 5\% percentile of activity bouth length
#' @return \code{percentile_activity_time_25} 25\% percentile of activity bout length
#' @return \code{percentile_activity_time_50} 50\% percentile of activity bout length
#' @return \code{percentile_activity_time_75} 75\% percentile of activity bout length
#' @return \code{percentile_activity_time_95} 95\% percentile of activity bout length
#' @return \code{alpha_activity} alpha of activity time, see details
#' @return \code{stepping_to_standing_ratio} Ratio of stepping to standing, the ratio of total stepping hours to standing hours
#' @details  \code{alpha_activity} is defined by \code{1+1/M}, where \code{M} is the average of \code{log(activity bout length /minimum activity bout length)}
#' @importFrom stats quantile
#' @examples
#' #For CRAN less than 5s running time policy, we only select the first day to run.
#' r4=SLI(sample_event[1:3095,],sample_bed_time[1,],sample_takeon_log[1,])
#' summary(r4)
#' @export
#'

SLI=function(final_dat,bed_time,takeoff_time) UseMethod("SLI")
#' @export

SLI.default=function(final_dat,bed_time,takeoff_time)
{
  out=SLI_summary(final_dat,bed_time,takeoff_time)
  out$call=match.call()
  class(out)="SLI"
  out
}
#' @export

print.SLI=function(x,...)
{cat("Call:\n")
  print(x$call)
}

#' @export

summary.SLI=function(object,...)
{
  TAB=cbind(object$year,object$month,object$day,object$dayofweek,object$weekday_or_weekend,object$total_number_of_activity_bouts,object$mean_activity_bout_length,object$prop_of_activity_time_greater_5min,object$prop_of_activity_time_greater_10min,object$prop_of_activity_time_greater_30min,object$total_activity_time_greater_5min,object$total_activity_time_greater_10min,object$total_activity_time_greater_30min,object$percentile_activity_time_5,object$percentile_activity_time_25,object$percentile_activity_time_50,object$percentile_activity_time_75,object$percentile_activity_time_95,object$alpha_activity,object$stepping_to_standing_ratio)

  colnames(TAB)=c("Year","Month","Day","Dayofweek","Weekday_or_weekend","total_number_of_activity_bouts","mean_activity_bout_length","prop_of_activity_time_greater_5min","prop_of_activity_time_greater_10min","prop_of_activity_time_greater_30min","total_activity_time_greater_5min","total_activity_time_greater_10min","total_activity_time_greater_30min","percentile_activity_time_5","percentile_activity_time_25","percentile_activity_time_50","percentile_activity_time_75","percentile_activity_time_95","alpha_activity","stepping_to_standing_ratio")
  res <- list(call=object$call,
              Table=TAB)
  class(res) <- "summary.SLI"
  res
}

