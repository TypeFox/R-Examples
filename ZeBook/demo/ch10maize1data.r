################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-14
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
library(maps)
# 1/ Exploration of available data on maize
# show data
maize.data_EuropeEU
list_n_sy=unique(maize.data_EuropeEU$sy)

GPS_sy=merge(maize.data_EuropeEU, unique(weather_EuropeEU[,c("idsite","GPSlatitude","GPSlongitude")]),by="idsite")

GPS_sy_B240=unique(GPS_sy[is.na(GPS_sy$LAIobs),c("idsite","GPSlatitude","GPSlongitude","year")])
GPS_sy_InSeasonBLAI=unique(GPS_sy[!is.na(GPS_sy$LAIobs),c("idsite","GPSlatitude","GPSlongitude","year")])
par(mar=c(0,0,0,0))
EuropeEU <- map(database = 'world', xlim=c(-11,28), ylim=c(34,63))
points(GPS_sy_B240$GPSlongitude,GPS_sy_B240$GPSlatitude, col = 'black', pch = 16,cex=1.5,)
points(GPS_sy_InSeasonBLAI$GPSlongitude,GPS_sy_InSeasonBLAI$GPSlatitude, col = 'darkgrey', pch = 16,cex=1.5)


# histogram of final biomass
hist(maize.data_EuropeEU[maize.data_EuropeEU$day==240,"Bobs"], xlab="Biomass observed at day=240", breaks=seq(2000,4000,by=250),main="")
# summary of Bobs at day=240
summary(maize.data_EuropeEU[maize.data_EuropeEU$day==240,"Bobs"])

# representation of dynamic variable for the k site-year
list_k_sy=unique(maize.data_EuropeEU[maize.data_EuropeEU$day!=240,"sy"])
plot(c(140,240),c(0,4000),type="n",xlab="day of year",ylab="observed biomass")
null=sapply(list_k_sy, function(sy) lines(maize.data_EuropeEU[maize.data_EuropeEU$sy==sy,c("day","Bobs")],type="b",lwd=2,pch=which(list_k_sy==sy)))

plot(c(140,240),c(0,8),type="n",xlab="day of year",ylab="observed LAI")
null=sapply(list_k_sy, function(sy) lines(maize.data_EuropeEU[maize.data_EuropeEU$sy==sy,c("day","LAIobs")],type="b",lwd=2,pch=which(list_k_sy==sy)))

# 2/ Exploration of available weather data
weather_all=weather_EuropeEU[(weather_EuropeEU$WEDAY>=100)&(weather_EuropeEU$WEDAY<=250),]
# function to compute sum of I and mean temperature between day 100 and day 250 (growing saison)
weather_analysis=function(sy,weather_all=weather_all) {
weather = maize.weather(working.year=strsplit(sy,"-")[[1]][2], working.site=strsplit(sy,"-")[[1]][1],weather_all=weather_all)

return(c(sumI=sum(weather$I), meanT=mean((weather$Tmin+weather$Tmax)/2)))}

# for the sample
wa_s=sapply(list_n_sy, weather_analysis, weather_all=weather_all)
# for all weather data site-year available
list_all_sy= apply(unique(weather_all[,c("idsite","WEYR")]),1,paste,collapse="-")
system.time(wa_all <- sapply(list_all_sy, weather_analysis, weather_all=weather_all))

# representation of all weather and the sample characteristics
h_sumI_all=hist(wa_all["sumI",], plot = FALSE, breaks=seq(1500,4000,by=500))
h_sumI_s=hist(wa_s["sumI",], plot = FALSE, breaks=seq(1500,4000,by=500))
h_meanT_all=hist(wa_all["meanT",],plot = FALSE, breaks=seq(10,30,by=5))
h_meanT_s=hist(wa_s["meanT",], plot = FALSE, breaks=seq(10,30,by=5))

h_sumI_all$counts=h_sumI_all$counts/sum(h_sumI_all$counts)
h_sumI_s$counts=h_sumI_s$counts/sum(h_sumI_s$counts)
h_meanT_all$counts=h_meanT_all$counts/sum(h_meanT_all$counts)
h_meanT_s$counts=h_meanT_s$counts/sum(h_meanT_s$counts)

par(mfrow=c(1,2))
plot( h_sumI_all, col="black", xlim=c(1500,4000), angle = -45,density = 20,main="",xlab="sum of radiation of growing season")
plot( h_sumI_s, col="black", xlim=c(1500,4000), add=T, angle = 45,density = 20)
legend("topright",c("all sy","sample"), angle = c(-45,45),density = 20,cex=0.70)
plot( h_meanT_all, col="black", xlim=c(10,30),angle = -45,density = 20,main="",xlab="mean temperature of growing season",ylim=c(0,0.8))
plot( h_meanT_s, col="black", xlim=c(10,30), add=T,angle = 45,density = 20)
legend("topright",c("all sy","sample"), angle = c(-45,45),density = 20,cex=0.70)

# end of file

