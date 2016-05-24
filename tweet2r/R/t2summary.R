#' @export
#' @import RSQLite ggplot2 plyr ggmap


t2summary<-function(tweets, geotweets){
  
  #------------------------------------------------------
  #Get crated_at as a a time stamp to work properly with it
  
  #change locale time language to avoid problems with posix
  
  #get local configuration and save
  lct<-Sys.getlocale("LC_TIME")
  #change language
  Sys.setlocale("LC_TIME", "C")
  
  #format createt_at column to timestamp
  list_time<-as.POSIXct(strptime(tweets$created_at, "%a %b %d %H:%M:%S %z %Y", tz = "UTC"), tz = "UTC")
 
  
  #back to previsous langue setting
  Sys.setlocale("LC_TIME", lct)
  #---------------------------------------------------------
  
  #---------------------------------------
  #create a data frame with the list_time
  
  #Separete the timestamp in year month day and hour
  year=strftime(list_time, "%Y")
  month=strftime(list_time, "%m")
  day=strftime(list_time, "%d")
  weekday=strftime(list_time, "%a")
  hour=strftime(list_time, "%H")
  data_time<-data.frame(year, month, day, hour)
  #-------------------------------------------------
  
  #-----------------------------------------------
  #count number of tweets per day and week day
  
  #----cont number of tweets per hour----------
  chour<-count(hour)
  #set column names
  colnames(chour)<-c("hour", "tweets")
  
  #plot tweets per hour using ggplot
  ghour<-ggplot(chour ,aes(hour,tweets))
  ghour<-ghour +geom_bar(stat="identity",fill="deepSkyblue2",colour = "deepSkyblue1")
  
  #cont number of tweets per weekday
  cweekday<-count(weekday)
  #set column names
  colnames(cweekday)<-c("day", "tweets")
  
  #plot tweets per weeday using ggplot
  gweekday<-ggplot(cweekday ,aes(day,tweets))
  gweekday <-gweekday +geom_bar(stat="identity",fill="deepSkyblue2",colour = "deepSkyblue1")
  
  #-----------------------------------------------
  #plot tweets on a map
  geotweets[,c("lon","lat")]
  
  #silly NULLS to pas chekcpackage
  lat=NULL
  lon=NULL
  mapt<-qmplot(lon,lat,data=geotweets@data, colour =I("deepSkyblue2"))
  
  #get number of tweets and geotweets
  ntweets<-nrow(tweets)
  ngeotweets<-nrow(geotweets)
  
  #get difference of tweets
  diftweets=ntweets-ngeotweets
  
  #get % of tweets geotweets
  pergeotweets=(ngeotweets/ntweets)*100
  
  summt<-data.frame(ntweets,ngeotweets,diftweets, pergeotweets)
  
  return(list("summt"=summt, "mapt"=mapt,"ghour"=ghour,"gweekday"=gweekday))

}


