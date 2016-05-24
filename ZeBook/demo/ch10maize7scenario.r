################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
# Utilization on scenarios
library(ZeBook)
library(maps)

list_all=unique(maize.weather(weather_all=weather_EuropeEU)[,c("idsite","GPSlatitude","GPSlongitude","year")])
list_all_sy=paste(list_all$idsite,list_all$year,sep="-")
list_all_site=unique(list_all[,c("idsite","GPSlatitude","GPSlongitude")])

# From estimation step, We keep RUE=2.175675
best_param <- maize.define.param()["nominal",]
param_opti=c("RUE")
best_param[param_opti]=2.175675

# prediction of mean B at 240 and std at 240.
system.time(sim<-maize.multisy(best_param, list_all_sy,sdate=100,ldate=250))
#save(sim,file="output/maize.case_05_sim.rda")
sim240=subset(sim, day==240)

matsy=matrix(as.numeric(unlist(strsplit(paste(sim240$sy),"-"))),ncol=2,byrow=TRUE)
sim240$idsite=matsy[,1]
sim240$year=matsy[,2]

mean_by_idsite<-as.matrix(unlist(by(sim240$B, sim240$idsite, mean)))
sd_by_idsite<-as.matrix(unlist(by(sim240$B, sim240$idsite, sd)))

by_idsite=data.frame(idsite=rownames(mean_by_idsite),B240.mean=mean_by_idsite ,B240.sd=sd_by_idsite )
by_idsite$B240.coefvar=by_idsite$B240.sd /  by_idsite$B240.mean
by_idsite=merge(list_all_site, by_idsite, by="idsite")

head(by_idsite)

################################################################################
# Mean B240
break.B240.mean= c(0 , 2250, 2750, 3000,3250,3750,+Inf)
class.m=cut(by_idsite$B240.mean, breaks=break.B240.mean, labels =c("<2250","2250-2750","2750-3000","3000-3250","3250-3750",">3750"))
col.m=paste(cut(by_idsite$B240.mean, breaks=break.B240.mean , labels =grey((seq(0.9,0.3,length=length(break.B240.mean)-1)))))

limits <- data.frame(limits = c('NO', 'NE', 'SO', 'SE'), lat = c(65, 65, 34, 34), long = c(-11, 30, -11, 30))
zone_select=list(x=c(-11,30,30,-11,-11), y=c(34,34,65,65,34))

# St/Cor B240
break.B240.coefvar= c(0 , 0.05, 0.1,0.25,0.5,1)
class.m.coefvar=cut(by_idsite$B240.coefvar, breaks=break.B240.coefvar, labels =c("<5%","5%-10%","10%-25%","25%-50%",">50%"))
col.m.coefvar=paste(cut(by_idsite$B240.coefvar, breaks=break.B240.coefvar , labels =grey((seq(0.9,0.3,length=length(break.B240.coefvar)-1)))))

limits <- data.frame(limits = c('NO', 'NE', 'SO', 'SE'), lat = c(65, 65, 34, 34), long = c(-11, 30, -11, 30))
zone_select=list(x=c(-11,30,30,-11,-11), y=c(34,34,65,65,34))

# Maps
# mean
par(mar=c(0.5,3,0.5,0.5))
EuropeEU <- map(database = 'world', xlim=c(limits$long[1],limits$long[2]), ylim=c(limits$lat[3],limits$lat[1]))
polygon(zone_select, col =rgb(0,0,0,0), border = rgb(0,0,0,0.2))
points(by_idsite$GPSlongitude,by_idsite$GPSlatitude, col= col.m,pch=16,cex=1)
legend("bottomright", legend=levels(class.m),  col=grey((seq(0.9,0.3,length=length(break.B240.mean)-1))), pch=rep(19,5),cex=0.8, bg="white")
# coefficient of variation
dev.new()
par(mar=c(0.5,3,0.5,0.5))
EuropeEU <- map(database = 'world', xlim=c(limits$long[1],limits$long[2]), ylim=c(limits$lat[3],limits$lat[1]))
polygon(zone_select, col =rgb(0,0,0,0), border = rgb(0,0,0,0.2))
points(by_idsite$GPSlongitude,by_idsite$GPSlatitude, col= col.m.coefvar,pch=16,cex=1)
legend("bottomright", legend=levels(class.m.coefvar),  col=grey((seq(0.9,0.3,length=length(break.B240.coefvar)-1))), pch=rep(19,5),cex=0.8, bg="white")

# end of file