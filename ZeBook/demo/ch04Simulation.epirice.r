################################################################################
# Working with dynamic models for agriculture
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2012-06-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
library(ZeBook)
# need to install the following package for the end of the demo.
#install.packages(c("maps","sp","spatial"))
library(maps)
library(sp)
library(spatial)

epirice.param=epirice.define.param()["nominal",]
weather_SouthAsia=weather_SouthAsia[,c("idsite","GPSlatitude","GPSlongitude","WEYR","WEDAY","TMAX","TMIN","RAIN","RH2M")]
################################################################################
# 1. Run model to find out the figure from Savary et al. (2012), used for validation.
# define a optimal weather for epidemy
weather=epirice.weather(working.year=1999, working.site=1)
weather$TMIN=10
weather$TMAX=30
weather$RH2M=95
weather$RAIN=0
# run model
simul=epirice.model(epirice.param, weather[90:300,],sdate=1,ldate=120,H0=600)
# compute AUDPC
AUDPC=sum(simul$severity,na.rm=TRUE)# %.day-1
AUDPC

# plot results
#png("epirice.optimal.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 12)
par(mar=c(3.1,3.1,1.1,3.1),cex=0.95, cex.axis=0.9)
#plot(simul$day,simul$TS,col="grey")
plot(range(simul$DACE),c(0,22000),type="n",xlab="",ylab="")
mtext("number of sites", side=2,line=2)
mtext("DACE (days after crop establishment)", side=1,line=2)
lines(simul$DACE,simul$H,col="black")
#lines(simul$DACE,simul$S,col="brown")
lines(simul$DACE,simul$L,col="orange",lwd=2,lty=2)
#lines(simul$DACE,simul$TOTDIS,col="darkred",lwd=2)
lines(simul$DACE,simul$II,col="red",lwd=2)
#lines(simul$day,simul$P,col="black")
legend("topleft",legend=c("H","L","II","severity"),lty=c(1,2,1,3),
lwd=c(1,2,2,3),col=c("black","orange","red","black"),
cex=0.6, ncol = 2)
par(ann=FALSE,new=TRUE)
plot.default(simul$DACE,simul$severity,axes=FALSE,type="l",xlim=range(simul$DACE),ylim=c(0,40),ylab='',lty=3, lwd=3)
axis(side=4)
mtext("severity (%)", side=4,line=2)
mtext(paste("AUDPC = ",round(AUDPC,1)),side=3,line=0,cex=0.8)
#dev.off()

################################################################################
# 2. Run model on multiple situations (sites and years) for a selected aera to obtain a map of disease (AUDPC)
date_transplantation = 150 # mean date of establishment for rice for North of India (June, 1th)
liste_site=unique(epirice.weather()[,c("idsite", "GPSlatitude","GPSlongitude")])
multi.simule= expand.grid(idsite=liste_site$id,year=unique(epirice.weather()$year), date_transplantation=date_transplantation)
# Warning : around 10 min of calculation
system.time(multi.simule <- epirice.multi.simule(epirice.param,multi.simule,all=TRUE))
multi.simule.mean_std=cbind(liste_site,AUDPC.mean=as.vector(by(multi.simule$AUDPC, multi.simule$idsite, mean)),AUDPC.sd=as.vector(by(multi.simule$AUDPC, multi.simule$idsite, sd)))
class.m=cut(multi.simule.mean_std$AUDPC.mean, breaks=c(0,2,8,25,57,361) , labels = c("1-2","3-8","9-25","26-57","58-361"))
col.m=paste(cut(multi.simule.mean_std$AUDPC.mean, breaks=c(0,2,8,25,57,361) , labels =grey((seq(0.9,0,length=5)))))

limits <- data.frame(limits = c('NO', 'NE', 'SO', 'SE'), lat = c(32, 32, 10, 10), long = c(78, 98, 78, 98))
dev.new()
#png("epirice_SouthAsia.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 12)
India <- map(database = 'world',c("India"),fill=TRUE, col="",type = "n", xlim=c(limits$long[1],limits$long[2]), ylim=c(limits$lat[3],limits$lat[1]), mar=c(0.2,0.2,0.2,0.2),  exact = TRUE)
polygon(India$x, India$y, col="white", border="black")
Bangladesh <- map(database = 'world',c("Bangladesh"),fill=TRUE, col="grey",type = "n",add=TRUE,  exact = TRUE)
polygon(Bangladesh$x, Bangladesh$y, col="white", border="black")
Myanmar <- map(database = 'world',c("Myanmar"),fill=TRUE, col="",type = "n",add=TRUE,  exact = TRUE)
polygon(Myanmar$x, Myanmar$y, col="white", border="black")
Nepal <- map(database = 'world',c("Nepal"),fill=TRUE, col="",type = "n",add=TRUE,  exact = TRUE)
polygon(Nepal$x, Nepal$y, col="white", border="black")

zone_select=list(x=c(min(liste_site$GPSlongitude),max(liste_site$GPSlongitude),
max(liste_site$GPSlongitude),min(liste_site$GPSlongitude),min(liste_site$GPSlongitude)),
y=c(min(liste_site$GPSlatitude),min(liste_site$GPSlatitude),max(liste_site$GPSlatitude),
max(liste_site$GPSlatitude),min(liste_site$GPSlatitude)))
polygon(zone_select, col =rgb(0,0,0,0), border = rgb(0,0,0,0.2))

# plot simple map
multi.simule.mean_std_t = multi.simule.mean_std[,c("GPSlongitude","GPSlatitude", "AUDPC.mean")]
points(multi.simule.mean_std$GPSlongitude,multi.simule.mean_std$GPSlatitude, col= col.m,pch=16,cex=1.5)
legend("bottomright", legend=unique(class.m),  col=unique(col.m), pch=rep(19,5))
#dev.off()

dev.new()
# Krigeage
India <- map(database = 'world',c("India"),fill=TRUE, col="",type = "n", xlim=c(limits$long[1],limits$long[2]), ylim=c(limits$lat[3],limits$lat[1]), mar=c(0.2,0.2,0.2,0.2),  exact = TRUE)
Bangladesh <- map(database = 'world',c("Bangladesh"),fill=TRUE, col="grey",type = "n",add=TRUE,  exact = TRUE)
Myanmar <- map(database = 'world',c("Myanmar"),fill=TRUE, col="",type = "n",add=TRUE,  exact = TRUE)
Nepal <- map(database = 'world',c("Nepal"),fill=TRUE, col="",type = "n",add=TRUE,  exact = TRUE)
#dev.off()

names(multi.simule.mean_std_t)=c("x","y","z")
multi.simule.mean_std.kr <- surf.ls(2, multi.simule.mean_std_t )
#variogram(multi.simule.mean_std.kr, 250)
trsurf <- trmat(multi.simule.mean_std.kr,  min(liste_site$GPSlongitude),max(liste_site$GPSlongitude) ,min(liste_site$GPSlatitude), max(liste_site$GPSlatitude), n=99)

trsurf$z[trsurf$z<0]=0

for (i in 1:dim(trsurf$z)[1]) {
  for (j in 1:dim(trsurf$z)[2]){
    pts=c(x=trsurf$x[i],y=trsurf$y[j])
    #print(pts)
    if ((((point.in.polygon(pts["x"],pts["y"], India$x, India$y, mode.checked=FALSE))|
(point.in.polygon(pts["x"],pts["y"], Bangladesh$x, Bangladesh$y, mode.checked=FALSE))|
(point.in.polygon(pts["x"],pts["y"], Myanmar$x, Myanmar$y, mode.checked=FALSE))|
(point.in.polygon(pts["x"],pts["y"], Nepal$x, Nepal$y, mode.checked=FALSE)))&
(point.in.polygon(pts["x"],pts["y"], zone_select$x, zone_select$y, mode.checked=FALSE)))==0)
    {
    #print("NA")
    trsurf$z[i,j]=NA
    }
  }}

# Combine Contour and map
#png("epirice_SouthAsia_krigeage.png", width = 8.5, height = 8.5,units = "cm", res = 300,pointsize = 6)
par(mar=c(0.5,2,0.5,0.5))
filled.contour(trsurf,  xlim=c(limits$long[1],limits$long[2]), ylim=c(limits$lat[3],limits$lat[1]),zlim=c(0,50), color = function(n){grey((seq(0.9,0,length=n)))},  levels=pretty(c(0,50), 10), asp=1,
plot.axes = {map(database = 'world',c("India","Nepal", "Bangladesh","Myanmar"),boundary = TRUE, add=TRUE,  exact = TRUE,  col="black")} )
#dev.off()

# end of file