#' Make hourly temperature record from daily data
#' 
#' This function generates hourly temperature records for a particular location
#' from daily minimum and maximum temperatures and latitude.
#' 
#' Temperature estimates are based on an idealized daily temperature curve that
#' uses a sine curve for daytime warming and a logarithmic decay function for
#' nighttime cooling. The input data frame can have more columns, which are
#' preserved, but ignored in the processing. References to papers outlining the
#' procedures are given below.
#' 
#' @param latitude the geographic latitude (in decimal degrees) of the location
#' of interest
#' @param year_file a data frame containing data on daily minimum temperature
#' (called Tmin), daily maximum temperature (called Tmax), and date
#' information. Dates can either be specified by two columns called Year and
#' JDay, which contain the Year and Julian date (day of the year), or as three
#' columns called Year, Month and Day. year_file cannot have any missing
#' values, so it may be a good idea to process the relevant columns with
#' make_all_day_table and interpolate_gaps before.
#' @return data frame containing all the columns of year_file, plus 24 columns
#' for hourly temperatures (called Hour_1 ... Hour_24).
#' @author Eike Luedeling
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' 
#' Luedeling E, Girvetz EH, Semenov MA and Brown PH, 2011. Climate change
#' affects winter chill for temperate fruit and nut trees. PLoS ONE 6(5),
#' e20155.
#' 
#' The temperature interpolation is described in
#' 
#' Linvill DE, 1990. Calculating chilling hours and chill units from daily
#' maximum and minimum temperature observations. HortScience 25(1), 14-16.
#' 
#' Calculation of sunrise, sunset and daylength was done according to
#' 
#' Spencer JW, 1971. Fourier series representation of the position of the Sun.
#' Search 2(5), 172.
#' 
#' Almorox J, Hontoria C and Benito M, 2005. Statistical validation of
#' daylength definitions for estimation of global solar radiation in Toledo,
#' Spain. Energy Conversion and Management 46(9-10), 1465-1471)
#' @keywords utility
#' @examples
#' 
#' weather<-fix_weather(KA_weather)
#' 
#' THourly<-make_hourly_temps(50.4,weather$weather)
#' 
#' #in most cases, you're probably better served by stack_hour_temperatures
#' 
#' @export make_hourly_temps
make_hourly_temps <-
function (latitude,year_file)

 {
         year_file<-year_file[which(!is.na(year_file$Tmin)&!is.na(year_file$Tmax)),]
  
         if(!"JDay" %in% colnames(year_file))
         year_file[,"JDay"]<-strptime(paste(year_file$Month,"/",year_file$Day,"/",year_file$Year,sep=""),"%m/%d/%Y")$yday+1
  
         preserve_columns<-colnames(year_file)
         year_file$Gamma<-2*pi/365*((year_file$JDay)-1)
         year_file$Delta<-180/pi*(0.006918-0.399912*cos(year_file$Gamma)+0.070257*sin(year_file$Gamma)-0.006758*cos(year_file$Gamma)+0.000907*sin(year_file$Gamma)-0.002697*cos(3*(year_file$Gamma))+0.00148*sin(3*(year_file$Gamma)))
         year_file$CosWo<-(sin(-0.8333/360*2*pi)-sin(latitude/360*2*pi)*sin(year_file$Delta/360*2*pi))/(cos(latitude/360*2*pi)*cos(year_file$Delta/360*2*pi))
         year_file$Sunrise[which(year_file$CosWo>=-1&year_file$CosWo<=1)]<-12-acos(year_file$CosWo)/(15/360*2*pi)
         if(length(which(year_file$CosWo<(-1)||year_file$CosWo>1))>0) {year_file$Sunrise[which(year_file$CosWo<-1||year_file$CosWo>1)]<--99}
         year_file$Sunset[which(year_file$CosWo>=-1&year_file$CosWo<=1)]<-12+acos(year_file$CosWo)/(15/360*2*pi)
         if(length(which(year_file$CosWo<(-1)||year_file$CosWo>1))>0) {year_file$Sunset[which(year_file$CosWo<-1||year_file$CosWo>1)]<--99}
         year_file$Daylength[which(year_file$CosWo>=-1&year_file$CosWo<=1)]<-2*acos(year_file$CosWo)/(15/360*2*pi)
         if(length(which(year_file$CosWo<(-1)||year_file$CosWo>1))>0) {year_file$Daylength[which(year_file$CosWo<-1||year_file$CosWo>1)]<--99}
         year_file$prev_max<-year_file$Tmax[c(nrow(year_file),1:(nrow(year_file)-1))]
         year_file$next_min<-year_file$Tmin[c(2:nrow(year_file),1)]
         year_file$prev_min<-year_file$Tmin[c(nrow(year_file),1:(nrow(year_file)-1))]
         year_file$Tsunset<-year_file$Tmin+(year_file$Tmax-year_file$Tmin)*
                            sin((pi*(year_file$Sunset-year_file$Sunrise)/(year_file$Daylength+4)))
         year_file$prev_Tsunset<-year_file$prev_min+(year_file$prev_max-year_file$prev_min)*
                            sin((pi*(year_file$Sunset-year_file$Sunrise)/(year_file$Daylength+4)))
         colnum<-ncol(year_file)+1

         hourcol<-c(colnum:(colnum+23))

     for (hourcount in 0:23)
             {

              if(length(which(year_file$Daylength==-99))>0) {year_file[which(year_file$Daylength==-99),colnum+hourcount]<-(year_file$Tmax+year_file$Tmin)/2}

              c_morn<-which(hourcount<=year_file$Sunrise)
              c_day<-which(hourcount>year_file$Sunrise&hourcount<year_file$Sunset+1)
              c_eve<-which(hourcount>=year_file$Sunset+1)
              nn<-colnum+hourcount

              if(length(which(year_file$Daylength>(-99)))>0)
                   {if(length(c_morn)>0)     #before sunrise
                          {year_file[c_morn,nn]<-
                                year_file$prev_Tsunset[c_morn]-  #prev temp at sunset
                               ((year_file$prev_Tsunset[c_morn]-year_file$Tmin[c_morn])/
                                  log(24-year_file$Daylength[c_morn])*
                                  log(hourcount+24-year_file$Sunset[c_morn]))}

                    if(length(c_day)>0)     #between sunrise and an hour after sunset
                          {year_file[c_day,colnum+hourcount]<-
                            year_file$Tmin[c_day]+
                            (year_file$Tmax[c_day]-year_file$Tmin[c_day])*
                            sin((pi*(hourcount-year_file$Sunrise[c_day])/
                                  (year_file$Daylength[c_day]+4)))}

                    if(length(c_eve)>0)                   #after sunset
                           {year_file[c_eve,colnum+hourcount]<-
                               year_file$Tsunset[c_eve]- #temp at sunset
                               ((year_file$Tsunset[c_eve]-year_file$next_min[c_eve])/
                                  log(24-year_file$Daylength[c_eve])*
                                  log(hourcount-year_file$Sunset[c_eve]))}

                 }}
                 colnames(year_file)[(ncol(year_file)-23):(ncol(year_file))]<-c(paste("Hour_",1:24,sep=""))
                 year_file<-year_file[,c(preserve_columns,paste("Hour_",1:24,sep=""))]
                 year_file[1,(ncol(year_file)-23):(ncol(year_file))][which(is.na(year_file[1,(ncol(year_file)-23):(ncol(year_file))]))]<-year_file[1,"Tmin"]
                 year_file[nrow(year_file),(ncol(year_file)-23):(ncol(year_file))][which(is.na(year_file[nrow(year_file),(ncol(year_file)-23):(ncol(year_file))]))]<-year_file[nrow(year_file),"Tmin"]
         
         return(year_file)
                 }
