

MVPA_summary=function (final_dat,bed_time,takeoff_time) {
  if(missing(takeoff_time)){
    final_dat=clean_time(final_dat)}else{
      final_dat=clean_time(final_dat,takeoff_time)
    }

  if(missing(bed_time)){
    mvpa_sporadic=NULL
    mvpa=NULL
    is_interval_valid=NULL


    temp_mat=final_dat



    temp_mat_for_activity<- temp_mat
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

    table1=c()
    ll=0
    ###################################################
    for (i in unique(combined_temp_mat_for_activity$new_date)){
      ll=ll+1
      temp_mat_oneday=final_dat[final_dat$new_date==i,]
      # temp_mat_oneday<-temp_mat_oneday[temp_mat_oneday[,1]>record.getup.time[ll] &temp_mat_oneday[,1]<record.sleep.time[ll], ]

      combined_temp_mat_for_activity_oneday=combined_temp_mat_for_activity[combined_temp_mat_for_activity$new_date==i,]
      #combined_temp_mat_for_activity_oneday<-combined_temp_mat_for_activity_oneday[combined_temp_mat_for_activity_oneday[,1]>record.getup.time[ll] &combined_temp_mat_for_activity_oneday[,1]<record.sleep.time[ll], ]

      time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))
      month<-as.numeric(format(time_char,"%m"))
      day<-as.numeric(format(time_char,"%d"))
      year<-as.numeric(format(time_char,"%Y"))
      dayofweek<-as.numeric(format(time_char,"%w"))
      weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)

      ###############
      ###############Calculation
      ###############
      ###
      temp_activity<-subset(combined_temp_mat_for_activity_oneday,combined_temp_mat_for_activity_oneday$ActivityCode==2)$Interval
      mvpa_sporadic_interval<- 1/4          ### 1 means 1 minute; 0_5 means 30 seconds; 1/6 means 10 seconds


      mvpa_1min_mat<-temp_mat_oneday
      start_time<-mvpa_1min_mat$date_time[1]
      end_time<-  mvpa_1min_mat$date_time[nrow(mvpa_1min_mat)]
      num_1min_interval<- trunc( (end_time-start_time)*24*60*60/(60*mvpa_sporadic_interval) )   ###point 1
      #(end_time-start_time) != sum(temp_mat$Interval)/60/60/24 if there is take off
      interval_length<-mvpa_sporadic_interval*60/60/60/24
      interval_1min_start<-start_time+interval_length*(1:(num_1min_interval))
      ##################
      mvpa_record_start_time<-c(start_time, interval_1min_start)
      mvpa_record_end_time<-c(interval_1min_start,end_time)
      ################## if there is take off, they won't be in combine_original_pseudo_mat
      combine_original_pseudo_mat<-do.call(rbind,sapply(1: length(mvpa_record_start_time),function(ll){
        temp_mat<-subset(mvpa_1min_mat,mvpa_1min_mat[,1]+mvpa_1min_mat[,2]/24/60/60>mvpa_record_start_time[ll] & mvpa_1min_mat[,1]<mvpa_record_end_time[ll] )
        if(nrow(temp_mat)==0) return (NULL)
        if(temp_mat[nrow(temp_mat),1]+(temp_mat[nrow(temp_mat),2]/24/60/60)> mvpa_record_end_time[ll])  ###if this activity is the last one and it surpass the take off log time
        {
          temp_mat[nrow(temp_mat),4]<-temp_mat[nrow(temp_mat),4]*(mvpa_record_end_time[ll]-temp_mat[nrow(temp_mat),1])/(temp_mat[nrow(temp_mat),2]/24/60/60)
          temp_mat[nrow(temp_mat),2]<- (mvpa_record_end_time[ll]-temp_mat[nrow(temp_mat),1])*24*60*60
        }
        if(temp_mat[1,1]<mvpa_record_start_time[ll])   ###if this activity is the first one and itis earlier than the take on log time
        {
          temp_mat[1,4]<-temp_mat[1,4]* (temp_mat[1,2]-(mvpa_record_start_time[ll]- temp_mat[1,1])*24*60*60)/temp_mat[1,2]
          temp_mat[1,2]<-temp_mat[1,2]-(mvpa_record_start_time[ll]- temp_mat[1,1])*24*60*60
          temp_mat[1,1]<-mvpa_record_start_time[ll]
        }

        return(cbind(temp_mat,ll))
      }, simplify = F)  )

      colnames(combine_original_pseudo_mat)<-c("date_time","Interval","ActivityCode", "METs","new_date","one_minute_interval")
      ############################################################
      ############################################################ step2 summary 1 min intervals
      ############################################################
      one_minute_collection<-by(combine_original_pseudo_mat,combine_original_pseudo_mat$one_minute_interval,function(s)c(min(s$date_time),sum(s$METs)*(60/mvpa_sporadic_interval),unique(s$one_minute_interval),sum(s$Interval),min(s$new_date)  )) ###point 3
      one_minute_mat<-do.call(rbind,one_minute_collection)
      one_minute_mat<-subset(one_minute_mat,one_minute_mat[,3]!=0 & one_minute_mat[,4]>(60*mvpa_sporadic_interval*0.9) & one_minute_mat[,4]<(60*mvpa_sporadic_interval*1.1)     )  ### one_minute_mat[,4] is the true length, it may not be exactly 30 second, can have a few seconds bias
      #########################
      if(trunc(nrow(one_minute_mat)/(10/mvpa_sporadic_interval))==0)  {
        ten_minute_vec<-rep(1,nrow(one_minute_mat))
        if(nrow(one_minute_mat)==1) ten_minute_mat<- data.frame(t(c(one_minute_mat[1:length(ten_minute_vec),],ten_minute_vec)))
        if(nrow(one_minute_mat)>1) ten_minute_mat<- data.frame(cbind(one_minute_mat[1:length(ten_minute_vec),],ten_minute_vec))
      } else
      {
        ten_minute_vec<-rep(1:trunc(nrow(one_minute_mat)/(10/mvpa_sporadic_interval)),each=  (10/mvpa_sporadic_interval)    )    ##### 30s to 10 min ###point 4
        ten_minute_mat<- data.frame(cbind(one_minute_mat[1:length(ten_minute_vec),],ten_minute_vec))
      }

      colnames(ten_minute_mat)<-c("date_time","mets","one_minute_interval","interval_length","new_date","ten_minute_interval")
      ############################################################
      ############################################################ step3 summary 10 min intervals
      ############################################################
      #### if in 10 minutes, 8 minutes have METs>3, it is MVPA bout; if less than 8minutes, they are counted as mvpa sporadic_
      is_mvpa<-function(s) if(s>=(8/mvpa_sporadic_interval)    ) return(1) else return(0)   ###point 5  ### this is for MVPA long bout
      ten_minute_collection<- data.frame(do.call(rbind,by(ten_minute_mat,ten_minute_mat$ten_minute_interval,function(s)c(min(s$date_time), is_mvpa(length(which(s$mets>=3))), mean(s$mets),length(which(s$mets>=3)), abs(max(s$date_time)-min(s$date_time)-sum(s$interval_length[1:(length(s$interval_length)-1)])/24/60/60 )   ))))
      colnames(ten_minute_collection)<-c("date_time","mvpa","mets","mvpa_sporadic","is_interval_valid")  #### is_interval_valid is to avoid the wear off during the day problem
      ten_minute_collection<-subset(ten_minute_collection,is_interval_valid<0.003) ###if the interval has 5 minutes take off, we do not take it

      ############################################################
      ############################################################ step4 MVPA information
      ############################################################
      #### total time
      Total_MVPA_Long_Bout_time<-nrow(subset(ten_minute_collection,mvpa==1))/6  ###by hours
      Total_mvpa_sporadic_time<-sum(subset(ten_minute_collection,mvpa_sporadic>0 & mvpa!=1)$mvpa_sporadic )/(60/mvpa_sporadic_interval) ###by hours ###point 6
      Total_MVPA_time<- Total_MVPA_Long_Bout_time+Total_mvpa_sporadic_time
      Total_light_time<- sum(temp_activity) /60/60-Total_MVPA_time
      #### Long Bouts+Sporadic_time runs
      Long_Bouts_and_Sporadic_run<- rle(ifelse( one_minute_mat[,2]>=3,1,0))
      Total_Number_of_MVPA_Long_Bouts_and_Sporadic<-  length(which(Long_Bouts_and_Sporadic_run$values==1))

      run_for_Long_Bouts_and_Sporadic_mvpa<- Long_Bouts_and_Sporadic_run$lengths[which(Long_Bouts_and_Sporadic_run$values==1)]/ (60/mvpa_sporadic_interval) ###by hours ###point 7

      if(Total_MVPA_time==0)  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2<-0 else  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2<-100*length(run_for_Long_Bouts_and_Sporadic_mvpa[run_for_Long_Bouts_and_Sporadic_mvpa>1/30])/Total_Number_of_MVPA_Long_Bouts_and_Sporadic
      if(Total_MVPA_time==0)  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5<-0 else  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5<-100*length(run_for_Long_Bouts_and_Sporadic_mvpa[run_for_Long_Bouts_and_Sporadic_mvpa>1/12])/Total_Number_of_MVPA_Long_Bouts_and_Sporadic
      if(Total_MVPA_time==0)  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10<-0 else Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10<-100*length(run_for_Long_Bouts_and_Sporadic_mvpa[run_for_Long_Bouts_and_Sporadic_mvpa>1/6])/Total_Number_of_MVPA_Long_Bouts_and_Sporadic


      #### Long Bouts Only runs
      run_for_mvpa<-rle(ten_minute_collection$mvpa)

      Total_Number_of_MVPA_Long_Bouts<-  length(which(run_for_mvpa$values==1))
      temp_mvpa_long_bout<-run_for_mvpa$lengths[which(run_for_mvpa$values==1)]/6 ###by hours
      if(Total_Number_of_MVPA_Long_Bouts==0) Mean_MVPA_Long_Bout_Length<-0 else Mean_MVPA_Long_Bout_Length<-mean(temp_mvpa_long_bout)
      ####
      if(Total_Number_of_MVPA_Long_Bouts==0)  Proportion_of_MVPA_Long_Bouts_greater_10<-0 else  Proportion_of_MVPA_Long_Bouts_greater_10<-100*length(temp_mvpa_long_bout[temp_mvpa_long_bout>1/6])/Total_Number_of_MVPA_Long_Bouts
      if(Total_Number_of_MVPA_Long_Bouts==0)  Proportion_of_MVPA_Long_Bouts_greater_20<-0 else  Proportion_of_MVPA_Long_Bouts_greater_20<-100*length(temp_mvpa_long_bout[temp_mvpa_long_bout>2/6])/Total_Number_of_MVPA_Long_Bouts


      ####################################################
      #################################################### MET_value
      ####################################################
      Highest_MET_value_15s<- max(one_minute_mat[,2])
      Highest_MET_value_10min<- max(ten_minute_collection[,3])
      ################################################ MET_value from MVPA
      Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa<- sum((one_minute_mat[,2]/60/60*one_minute_mat[,4])[one_minute_mat[,2]>=3])
      Total_MET_hrs_Long_Bouts<- sum((subset(ten_minute_collection,mvpa==1))$mets/60*10)
      table<- cbind(year,month,day,dayofweek,weekday_or_weekend,Total_light_time,Total_MVPA_time, Total_MVPA_Long_Bout_time,Total_mvpa_sporadic_time,Total_Number_of_MVPA_Long_Bouts_and_Sporadic,Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2,Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5,Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10, Total_Number_of_MVPA_Long_Bouts,Mean_MVPA_Long_Bout_Length, Proportion_of_MVPA_Long_Bouts_greater_10,Proportion_of_MVPA_Long_Bouts_greater_20, Highest_MET_value_15s,Highest_MET_value_10min, Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa,Total_MET_hrs_Long_Bouts)

      row.names(table)=NULL
      table1=rbind(table1,table)
      #print(table1)
    }
  }else{
    mvpa_sporadic=NULL
    mvpa=NULL
    is_interval_valid=NULL


    temp_mat=final_dat



    temp_mat_for_activity<- temp_mat
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

    table1=c()
    ll=0
    ###################################################
    for (i in unique(combined_temp_mat_for_activity$new_date)){
      ll=ll+1
      temp_mat_oneday=final_dat[final_dat$new_date==i,]
      temp_mat_oneday<-temp_mat_oneday[temp_mat_oneday[,1]>record.getup.time[ll] &temp_mat_oneday[,1]<record.sleep.time[ll], ]

      combined_temp_mat_for_activity_oneday=combined_temp_mat_for_activity[combined_temp_mat_for_activity$new_date==i,]
      combined_temp_mat_for_activity_oneday<-combined_temp_mat_for_activity_oneday[combined_temp_mat_for_activity_oneday[,1]>record.getup.time[ll] &combined_temp_mat_for_activity_oneday[,1]<record.sleep.time[ll], ]

      time_char<-as.POSIXlt(i*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0))
      month<-as.numeric(format(time_char,"%m"))
      day<-as.numeric(format(time_char,"%d"))
      year<-as.numeric(format(time_char,"%Y"))
      dayofweek<-as.numeric(format(time_char,"%w"))
      weekday_or_weekend<- ifelse( dayofweek!=0 & dayofweek!=6,1,0)

      ###############
      ###############Calculation
      ###############
      ###
      temp_activity<-subset(combined_temp_mat_for_activity_oneday,combined_temp_mat_for_activity_oneday$ActivityCode==2)$Interval
      mvpa_sporadic_interval<- 1/4          ### 1 means 1 minute; 0_5 means 30 seconds; 1/6 means 10 seconds


      mvpa_1min_mat<-temp_mat_oneday
      start_time<-mvpa_1min_mat$date_time[1]
      end_time<-  mvpa_1min_mat$date_time[nrow(mvpa_1min_mat)]
      num_1min_interval<- trunc( (end_time-start_time)*24*60*60/(60*mvpa_sporadic_interval) )   ###point 1
      #(end_time-start_time) != sum(temp_mat$Interval)/60/60/24 if there is take off
      interval_length<-mvpa_sporadic_interval*60/60/60/24
      interval_1min_start<-start_time+interval_length*(1:(num_1min_interval))
      ##################
      mvpa_record_start_time<-c(start_time, interval_1min_start)
      mvpa_record_end_time<-c(interval_1min_start,end_time)
      ################## if there is take off, they won't be in combine_original_pseudo_mat
      combine_original_pseudo_mat<-do.call(rbind,sapply(1: length(mvpa_record_start_time),function(ll){
        temp_mat<-subset(mvpa_1min_mat,mvpa_1min_mat[,1]+mvpa_1min_mat[,2]/24/60/60>mvpa_record_start_time[ll] & mvpa_1min_mat[,1]<mvpa_record_end_time[ll] )
        if(nrow(temp_mat)==0) return (NULL)
        if(temp_mat[nrow(temp_mat),1]+(temp_mat[nrow(temp_mat),2]/24/60/60)> mvpa_record_end_time[ll])  ###if this activity is the last one and it surpass the take off log time
        {
          temp_mat[nrow(temp_mat),4]<-temp_mat[nrow(temp_mat),4]*(mvpa_record_end_time[ll]-temp_mat[nrow(temp_mat),1])/(temp_mat[nrow(temp_mat),2]/24/60/60)
          temp_mat[nrow(temp_mat),2]<- (mvpa_record_end_time[ll]-temp_mat[nrow(temp_mat),1])*24*60*60
        }
        if(temp_mat[1,1]<mvpa_record_start_time[ll])   ###if this activity is the first one and itis earlier than the take on log time
        {
          temp_mat[1,4]<-temp_mat[1,4]* (temp_mat[1,2]-(mvpa_record_start_time[ll]- temp_mat[1,1])*24*60*60)/temp_mat[1,2]
          temp_mat[1,2]<-temp_mat[1,2]-(mvpa_record_start_time[ll]- temp_mat[1,1])*24*60*60
          temp_mat[1,1]<-mvpa_record_start_time[ll]
        }

        return(cbind(temp_mat,ll))
      }, simplify = F)  )

      colnames(combine_original_pseudo_mat)<-c("date_time","Interval","ActivityCode", "METs","new_date","one_minute_interval")
      ############################################################
      ############################################################ step2 summary 1 min intervals
      ############################################################
      one_minute_collection<-by(combine_original_pseudo_mat,combine_original_pseudo_mat$one_minute_interval,function(s)c(min(s$date_time),sum(s$METs)*(60/mvpa_sporadic_interval),unique(s$one_minute_interval),sum(s$Interval),min(s$new_date)  )) ###point 3
      one_minute_mat<-do.call(rbind,one_minute_collection)
      one_minute_mat<-subset(one_minute_mat,one_minute_mat[,3]!=0 & one_minute_mat[,4]>(60*mvpa_sporadic_interval*0.9) & one_minute_mat[,4]<(60*mvpa_sporadic_interval*1.1)     )  ### one_minute_mat[,4] is the true length, it may not be exactly 30 second, can have a few seconds bias
      #########################
      if(trunc(nrow(one_minute_mat)/(10/mvpa_sporadic_interval))==0)  {
        ten_minute_vec<-rep(1,nrow(one_minute_mat))
        if(nrow(one_minute_mat)==1) ten_minute_mat<- data.frame(t(c(one_minute_mat[1:length(ten_minute_vec),],ten_minute_vec)))
        if(nrow(one_minute_mat)>1) ten_minute_mat<- data.frame(cbind(one_minute_mat[1:length(ten_minute_vec),],ten_minute_vec))
      } else
      {
        ten_minute_vec<-rep(1:trunc(nrow(one_minute_mat)/(10/mvpa_sporadic_interval)),each=  (10/mvpa_sporadic_interval)    )    ##### 30s to 10 min ###point 4
        ten_minute_mat<- data.frame(cbind(one_minute_mat[1:length(ten_minute_vec),],ten_minute_vec))
      }

      colnames(ten_minute_mat)<-c("date_time","mets","one_minute_interval","interval_length","new_date","ten_minute_interval")
      ############################################################
      ############################################################ step3 summary 10 min intervals
      ############################################################
      #### if in 10 minutes, 8 minutes have METs>3, it is MVPA bout; if less than 8minutes, they are counted as mvpa sporadic_
      is_mvpa<-function(s) if(s>=(8/mvpa_sporadic_interval)    ) return(1) else return(0)   ###point 5  ### this is for MVPA long bout
      ten_minute_collection<- data.frame(do.call(rbind,by(ten_minute_mat,ten_minute_mat$ten_minute_interval,function(s)c(min(s$date_time), is_mvpa(length(which(s$mets>=3))), mean(s$mets),length(which(s$mets>=3)), abs(max(s$date_time)-min(s$date_time)-sum(s$interval_length[1:(length(s$interval_length)-1)])/24/60/60 )   ))))
      colnames(ten_minute_collection)<-c("date_time","mvpa","mets","mvpa_sporadic","is_interval_valid")  #### is_interval_valid is to avoid the wear off during the day problem
      ten_minute_collection<-subset(ten_minute_collection,is_interval_valid<0.003) ###if the interval has 5 minutes take off, we do not take it

      ############################################################
      ############################################################ step4 MVPA information
      ############################################################
      #### total time
      Total_MVPA_Long_Bout_time<-nrow(subset(ten_minute_collection,mvpa==1))/6  ###by hours
      Total_mvpa_sporadic_time<-sum(subset(ten_minute_collection,mvpa_sporadic>0 & mvpa!=1)$mvpa_sporadic )/(60/mvpa_sporadic_interval) ###by hours ###point 6
      Total_MVPA_time<- Total_MVPA_Long_Bout_time+Total_mvpa_sporadic_time
      Total_light_time<- sum(temp_activity) /60/60-Total_MVPA_time
      #### Long Bouts+Sporadic_time runs
      Long_Bouts_and_Sporadic_run<- rle(ifelse( one_minute_mat[,2]>=3,1,0))
      Total_Number_of_MVPA_Long_Bouts_and_Sporadic<-  length(which(Long_Bouts_and_Sporadic_run$values==1))

      run_for_Long_Bouts_and_Sporadic_mvpa<- Long_Bouts_and_Sporadic_run$lengths[which(Long_Bouts_and_Sporadic_run$values==1)]/ (60/mvpa_sporadic_interval) ###by hours ###point 7

      if(Total_MVPA_time==0)  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2<-0 else  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2<-100*length(run_for_Long_Bouts_and_Sporadic_mvpa[run_for_Long_Bouts_and_Sporadic_mvpa>1/30])/Total_Number_of_MVPA_Long_Bouts_and_Sporadic
      if(Total_MVPA_time==0)  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5<-0 else  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5<-100*length(run_for_Long_Bouts_and_Sporadic_mvpa[run_for_Long_Bouts_and_Sporadic_mvpa>1/12])/Total_Number_of_MVPA_Long_Bouts_and_Sporadic
      if(Total_MVPA_time==0)  Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10<-0 else Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10<-100*length(run_for_Long_Bouts_and_Sporadic_mvpa[run_for_Long_Bouts_and_Sporadic_mvpa>1/6])/Total_Number_of_MVPA_Long_Bouts_and_Sporadic


      #### Long Bouts Only runs
      run_for_mvpa<-rle(ten_minute_collection$mvpa)

      Total_Number_of_MVPA_Long_Bouts<-  length(which(run_for_mvpa$values==1))
      temp_mvpa_long_bout<-run_for_mvpa$lengths[which(run_for_mvpa$values==1)]/6 ###by hours
      if(Total_Number_of_MVPA_Long_Bouts==0) Mean_MVPA_Long_Bout_Length<-0 else Mean_MVPA_Long_Bout_Length<-mean(temp_mvpa_long_bout)
      ####
      if(Total_Number_of_MVPA_Long_Bouts==0)  Proportion_of_MVPA_Long_Bouts_greater_10<-0 else  Proportion_of_MVPA_Long_Bouts_greater_10<-100*length(temp_mvpa_long_bout[temp_mvpa_long_bout>1/6])/Total_Number_of_MVPA_Long_Bouts
      if(Total_Number_of_MVPA_Long_Bouts==0)  Proportion_of_MVPA_Long_Bouts_greater_20<-0 else  Proportion_of_MVPA_Long_Bouts_greater_20<-100*length(temp_mvpa_long_bout[temp_mvpa_long_bout>2/6])/Total_Number_of_MVPA_Long_Bouts


      ####################################################
      #################################################### MET_value
      ####################################################
      Highest_MET_value_15s<- max(one_minute_mat[,2])
      Highest_MET_value_10min<- max(ten_minute_collection[,3])
      ################################################ MET_value from MVPA
      Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa<- sum((one_minute_mat[,2]/60/60*one_minute_mat[,4])[one_minute_mat[,2]>=3])
      Total_MET_hrs_Long_Bouts<- sum((subset(ten_minute_collection,mvpa==1))$mets/60*10)
      table<- cbind(year,month,day,dayofweek,weekday_or_weekend,Total_light_time,Total_MVPA_time, Total_MVPA_Long_Bout_time,Total_mvpa_sporadic_time,Total_Number_of_MVPA_Long_Bouts_and_Sporadic,Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2,Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5,Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10, Total_Number_of_MVPA_Long_Bouts,Mean_MVPA_Long_Bout_Length, Proportion_of_MVPA_Long_Bouts_greater_10,Proportion_of_MVPA_Long_Bouts_greater_20, Highest_MET_value_15s,Highest_MET_value_10min, Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa,Total_MET_hrs_Long_Bouts)

      row.names(table)=NULL
      table1=rbind(table1,table)
      #print(table1)
    }
  }
  #colnames(table1)<- c("year","month","day","dayofweek","weekday_or_weekend","Total_light_time","Total_MVPA_time","Total_MVPA_Long_Bout_time","Total_mvpa_sporadic_time","Total_Number_of_MVPA_Long_Bouts_and_Sporadic","Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2","Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5","Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10","Total_Number_of_MVPA_Long_Bouts","Mean_MVPA_Long_Bout_Length","Proportion_of_MVPA_Long_Bouts_greater_10","Proportion_of_MVPA_Long_Bouts_greater_20","Highest_MET_value_15s","Highest_MET_value_10min","Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa","Total_MET_hrs_Long_Bouts")

  #out=list( total_number_of_activity_bouts=total_number_of_activity_bouts,mean_activity_bout_length=mean_activity_bout_length,prop_of_activity_time_greater_5min=prop_of_activity_time_greater_5min,prop_of_activity_time_greater_10min=prop_of_activity_time_greater_10min,prop_of_activity_time_greater_30min=prop_of_activity_time_greater_30min,total_activity_time_greater_5min=total_activity_time_greater_5min,total_activity_time_greater_10min=total_activity_time_greater_10min,total_activity_time_greater_30min=total_activity_time_greater_30min,percentile_activity_time_5=percentile_activity_time_5,percentile_activity_time_25=percentile_activity_time_25,percentile_activity_time_50=percentile_activity_time_50,percentile_activity_time_75=percentile_activity_time_75,percentile_activity_time_95=percentile_activity_time_95,alpha_activity=alpha_activity,stepping_to_standing_ratio=stepping_to_standing_ratio,table=table)
  #return(table1)

  list(year=table1[,1],month=table1[,2], day=table1[,3],dayofweek=table1[,4], weekday_or_weekend=table1[,5],Total_light_time=table1[,6], Total_MVPA_time=table1[,7], Total_MVPA_Long_Bout_time=table1[,8], Total_mvpa_sporadic_time=table1[,9],Total_Number_of_MVPA_Long_Bouts_and_Sporadic=table1[,10],Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2=table1[,11],Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5=table1[,12],Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10=table1[,13],Total_Number_of_MVPA_Long_Bouts=table1[,14],Mean_MVPA_Long_Bout_Length=table1[,15],Proportion_of_MVPA_Long_Bouts_greater_10=table1[,16],Proportion_of_MVPA_Long_Bouts_greater_20=table1[,17],Highest_MET_value_15s=table1[,18],Highest_MET_value_10min=table1[,19],Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa=table1[,20],Total_MET_hrs_Long_Bouts=table1[,21])

}
#' Moderate to Vigorous Physical Activity Summary
#'
#' Summarize moderate to vigorous activity measures
#' @param final_dat Raw event file, will be cleaned in this function. Event file is required for this function.
#' @param takeoff_time Take on and off time log, reported by participants. Log is not required for this function.
#' @param bed_time Sleep and wake up time log, reported by participants. Log is not required for this function.
#' @return \code{Year} The calendar year of recorded event
#' @return \code{Month} The calendar month of recorded event
#' @return \code{Day}   The calendar day of recorded event
#' @return \code{Dayofweek} The day of that week
#' @return \code{Weekday_or_weekend} The recored event date is a weekday or weekend (0 for weekday, 1 for weekend)
#' @return \code{Total_light_time} Total light time: total activity minutes minus the total MVPA hours
#' @return \code{Total_MVPA_time}  Total MVPA hours: summation of the total long bouts MVPA and total sporadic MVPA durations
#' @return \code{Total_MVPA_Long_Bout_time} Total MVPA long bout time: summation of the total long bouts MVPA durations
#' @return \code{Total_mvpa_sporadic_time} Total MVPA sporadic time: summation of the total sporadic MVPA durations
#' @return \code{Total_Number_of_MVPA_Long_Bouts_and_Sporadic} Total Number of MVPA (long bouts or sporadic): count the number of MVPA from 15 second data
#' @return \code{Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2} Proportion of MVPA (long bouts or sporadic) greater than 2 minutes
#' @return \code{Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5} Proportion of MVPA (long bouts or sporadic) greater than 5 minutes
#' @return \code{Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10} Proportion of MVPA (long bouts or sporadic) greater than 10 minutes
#' @return \code{Total_Number_of_MVPA_Long_Bouts} Total number of long bouts MVPA: count the number of long bouts MVPA from 10 minute data
#' @return \code{Mean_MVPA_Long_Bout_Length} Mean of long bout MVPA length
#' @return \code{Proportion_of_MVPA_Long_Bouts_greater_10} Proportion of long bouts MVPA greater than 10 minutes
#' @return \code{Proportion_of_MVPA_Long_Bouts_greater_20} Proportion of long bouts MVPA greater than 20 minutes
#' @return \code{Highest_MET_value_15s} Highest METs values in 15 second
#' @return \code{Highest_MET_value_10min} Highest METs values in 10 minutes
#' @return \code{Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa} Total METs hours from long bouts MVPA and sporadic MVPA: the summation of METs hours from all MVPA records
#' @return \code{Total_MET_hrs_Long_Bouts} Total METs hours from long bouts MVPA: the summation of METs hours from long bouts MVPA records
#' @details MVPA is defined into two types: long bout MVPA and sporadic MVPA. Long bout MVPA is defined as 10 consecutive minutes with METs>=3 (allowing 2 min below that threshold). Sporadic MVPA is defined as activities at any time with METS>=3 and they are not in long bouts MVPA.
#' @details Highest METs values in 15 second/10 minutes are calculated by picking up the maximum METs values from the combined data with 15 second intervals and the data with 10 minutes intervals, respectively.
#' @details Function \code{MVPA()} takes longer running time than other functions as its nature of more complex calculation.
#' @examples
#' #For CRAN less than 5s running time policy, we only select the first day to run.
#' r3=MVPA(sample_event[1:3095,],sample_bed_time[1,],sample_takeon_log[1,])
#' summary(r3)
#' @export
MVPA=function(final_dat,bed_time,takeoff_time) UseMethod("MVPA")

#' @export

MVPA.default=function(final_dat,bed_time,takeoff_time)
{
  out=MVPA_summary(final_dat,bed_time,takeoff_time)
  out$call=match.call()
  class(out)="MVPA"
  out
}

#' @export

print.MVPA=function(x,...)
{cat("Call:\n")
  print(x$call)
}

#' @export

summary.MVPA=function(object,...)
{
  TAB=cbind(object$year,object$month,object$day,object$dayofweek,object$weekday_or_weekend,object$Total_light_time,object$Total_MVPA_time,object$Total_MVPA_Long_Bout_time,object$Total_mvpa_sporadic_time,object$Total_Number_of_MVPA_Long_Bouts_and_Sporadic,object$Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2,object$Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5,object$Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10,object$Total_Number_of_MVPA_Long_Bouts,object$Mean_MVPA_Long_Bout_Length,object$Proportion_of_MVPA_Long_Bouts_greater_10,object$Proportion_of_MVPA_Long_Bouts_greater_20,object$Highest_MET_value_15s,object$Highest_MET_value_10min,object$Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa,object$Total_MET_hrs_Long_Bouts)

  colnames(TAB)=c("Year","Month","Day","Dayofweek","Weekday_or_weekend","Total_light_time","Total_MVPA_time","Total_MVPA_Long_Bout_time","Total_mvpa_sporadic_time","Total_Number_of_MVPA_Long_Bouts_and_Sporadic","Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_2","Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_5","Proportion_of_MVPA_Long_Bouts_and_Sporadic_greater_10","Total_Number_of_MVPA_Long_Bouts","Mean_MVPA_Long_Bout_Length","Proportion_of_MVPA_Long_Bouts_greater_10","Proportion_of_MVPA_Long_Bouts_greater_20","Highest_MET_value_15s","Highest_MET_value_10min","Total_MET_hrs_Long_Bouts_and_Sporadic_mvpa","Total_MET_hrs_Long_Bouts")
  res <- list(call=object$call,
              Table=TAB)
  class(res) <- "summary.MVPA"
  res
}


