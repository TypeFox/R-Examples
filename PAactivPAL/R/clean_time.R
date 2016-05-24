clean_time=function(dat,takeoff_time){
  if(missing(takeoff_time))
  {if(is.numeric(dat$Time)==F) ####if  is.numeric(dat$Time)==F, we need further modification of time in next step
  {
    ee<-as.character(dat$Time)
    max_length<-max(nchar(ee))
    ee[nchar(ee)!=max_length]<-"#1899-12-30 00:00:00#"
    #### use character type, may not be good
    ee_new<- (as.numeric( as.POSIXlt( substr(ee, 2, max_length-1 )  )    )+2209136400)/24/60/60
    #### use interval type, this is the best
    start_ee<-  min(which(dat[,2]>0))-1
    if (start_ee>1)
    {
      ee_new_int_type <-c( ee_new[1:((start_ee)-1)],ee_new[start_ee]+(dat$DataCount[start_ee:nrow(dat)]/10/24/60/60)  )
    }else { start_ee=1
      ee_new_int_type <-ee_new[start_ee]+(dat$DataCount[start_ee:nrow(dat)]/10/24/60/60)
    } #### if interval type has large difference with character type, use character type
    int_dif_char<-which(abs(ee_new-ee_new_int_type)>0.1 )
    ee_new_int_type[int_dif_char]<-ee_new[int_dif_char]
    ####
    dat<-cbind(ee_new_int_type,dat[,2:6])
  }
    final_dat<-dat[,c(1,3,4,6)]
    colnames(final_dat)<-c("date_time","Interval","ActivityCode", "METs")
    final_dat=final_dat[final_dat$date_time!=0,] #delete rows with time as #1899-12-30 00:00:00#

    new_date=0
    final_dat=cbind(final_dat,new_date)
    for (i in 1:nrow(final_dat)){
      #eg: 2010-03-31 23:59:59 is 40268.9583217593; 2010-04-01 00:00:00 is 40268.9583333333 but this is count as 2010-03-31
      if(final_dat$date_time[i]-floor(final_dat$date_time[i])<=0.125){ #eg: 2010-04-01 04:00:00 is 40269.125 but it is counted as 2010-03-31 so overwrite the date using 2010-03-31
        final_dat$new_date[i]= floor(final_dat$date_time[i])-1+0.125 }else{
          final_dat$new_date[i]= floor(final_dat$date_time[i])+0.125}
    }


    return(final_dat)}else{
  record<-takeoff_time

  ##
 if(nrow(record)==0)next
  ##
  record.start.time<-c()
  record.end.time<-c()
  for (kk in 1:nrow(record))
  {
    temp.start.time<-as.numeric(as.POSIXlt(strptime(as.character(paste(as.character(record[kk,3] ),as.character(record[kk,4]))),"%m/%d/%Y %H:%M:%S"))+2209136400)/24/60/60
    record.start.time<-c(record.start.time,temp.start.time)

    temp.end.time<-as.numeric(as.POSIXlt( strptime(as.character(paste(as.character(record[kk,5] ),as.character(record[kk,6]))),"%m/%d/%Y %H:%M:%S")  )+2209136400)/24/60/60
    record.end.time<-c(record.end.time,temp.end.time)
  }

  ##########################################  match get up time with take on time/// if multiple take on time in a day, the second take on time is seen as get up time
  #####list days take on
  #####time.char<-as.numeric(format(as.POSIXlt(record.start.time*24*60*60, origin = ISOdatetime(1899,12,30,0,0,0) ),"%d"))
  ###
#   record.getup<-subset(bed.time,bed.time$id==person)
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
if(is.numeric(dat$Time)==F) ####if  is.numeric(dat$Time)==F, we need further modification of time in next step
{
  ee<-as.character(dat$Time)
  max_length<-max(nchar(ee))
  ee[nchar(ee)!=max_length]<-"#1899-12-30 00:00:00#"
  #### use character type, may not be good
  ee_new<- (as.numeric( as.POSIXlt( substr(ee, 2, max_length-1 )  )    )+2209136400)/24/60/60
  #### use interval type, this is the best
  start_ee<-  min(which(dat[,2]>0))-1
  if (start_ee>1)
  {
    ee_new_int_type <-c( ee_new[1:((start_ee)-1)],ee_new[start_ee]+(dat$DataCount[start_ee:nrow(dat)]/10/24/60/60)  )
  }else
  {    start_ee=1
  ee_new_int_type <-ee_new[start_ee]+(dat$DataCount[start_ee:nrow(dat)]/10/24/60/60)
  } #### if interval type has large difference with character type, use character type
  int_dif_char<-which(abs(ee_new-ee_new_int_type)>0.1 )
  ee_new_int_type[int_dif_char]<-ee_new[int_dif_char]
  ####
  dat<-cbind(ee_new_int_type,dat[,2:6])
}
final_dat<-dat[,c(1,3,4,6)]
colnames(final_dat)<-c("date_time","Interval","ActivityCode", "METs")
final_dat=final_dat[final_dat$date_time!=0,] #delete rows with time as #1899-12-30 00:00:00#
if(!is.null(takeoff_time)){
final_dat<-do.call(rbind,sapply( 1: length(record.start.time),function(ll){
temp.mat<-final_dat[(final_dat[,1]+final_dat[,2]/24/60/60)>record.start.time[ll] & final_dat[,1]<record.end.time[ll] ,]
  if(nrow(temp.mat)==0) return (NULL)
  if(temp.mat[nrow(temp.mat),1]+(temp.mat[nrow(temp.mat),2]/24/60/60)> record.end.time[ll])  ###if this activity is the last one and it surpass the take off log time
  {
    temp.mat[nrow(temp.mat),4]<-temp.mat[nrow(temp.mat),4]*(record.end.time[ll]-temp.mat[nrow(temp.mat),1])/(temp.mat[nrow(temp.mat),2]/24/60/60)
    temp.mat[nrow(temp.mat),2]<- (record.end.time[ll]-temp.mat[nrow(temp.mat),1])*24*60*60
  }
  if(temp.mat[1,1]<record.start.time[ll])   ###if this activity is the first one and itis earlier than the take on log time
  {
    temp.mat[1,4]<-temp.mat[1,4]* (temp.mat[1,2]-(record.start.time[ll]- temp.mat[1,1])*24*60*60)/temp.mat[1,2]
    temp.mat[1,2]<-temp.mat[1,2]-(record.start.time[ll]- temp.mat[1,1])*24*60*60
    temp.mat[1,1]<-record.start.time[ll]
  }

  return(temp.mat)
}, simplify = F))
}
new_date=0
final_dat=cbind(final_dat,new_date)
for (i in 1:nrow(final_dat)){
  #eg: 2010-03-31 23:59:59 is 40268.9583217593; 2010-04-01 00:00:00 is 40268.9583333333 but this is count as 2010-03-31
  if(final_dat$date_time[i]-floor(final_dat$date_time[i])<=0.125){ #eg: 2010-04-01 04:00:00 is 40269.125 but it is counted as 2010-03-31 so overwrite the date using 2010-03-31
    final_dat$new_date[i]= floor(final_dat$date_time[i])-1+0.125 }else{
      final_dat$new_date[i]= floor(final_dat$date_time[i])+0.125}
}


return(final_dat)
}
}
