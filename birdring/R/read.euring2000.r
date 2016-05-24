# read EURING data format EURING 2000 into R
# Author: Fraenzi Korner-Nievergelt
#
# Reference: Speek et al. 2007:The EURING exchange-code 2000. www.euring.org
#-----------------------------------------------------------------------
# report of changes
# 22. 1. 2014 fk: condition and circumstances were transformed into numeric variables

#-------------------------------------------------------------------------------
read.euring2000 <- function(filename){
#-------------------------------------------------------------------------------
# filename: name of the txt-file obtained from EURING
#------------------------------------------------------------------------------- 
options(encoding="latin1")
warnorig <- options("warn")
options(warn= -1)

rodat<-read.table(filename, colClasses = "character", sep=",")
dat<-data.frame(scheme=substr(rodat$V1, 1,3), id.method=substr(rodat$V1, 4, 5), ring=substr(rodat$V1, 6, 15))
dat$ring.verif<-substr(rodat$V1, 16, 16)
dat$metal.ring.info<-substr(rodat$V1, 17, 17)
dat$marks.info<-substr(rodat$V1, 18, 19)
dat$spec.byringer<-substr(rodat$V1, 20, 24)
dat$spec.byscheme<-substr(rodat$V1, 25, 29)
dat$manipulated<-substr(rodat$V1, 30, 30)
dat$moved<-substr(rodat$V1, 31, 31)
dat$catching.method<-substr(rodat$V1, 32, 32)
dat$catching.lures<-substr(rodat$V1, 33, 33)
dat$sex.byringer<-substr(rodat$V1, 34, 34)
dat$sex.byscheme<-substr(rodat$V1, 35, 35)
dat$age.byringer<-substr(rodat$V1, 36, 36)
dat$age.byscheme<-substr(rodat$V1, 37, 37)
dat$status<-substr(rodat$V1, 38, 38)
dat$broodsize<-substr(rodat$V1, 39, 40)
dat$pullus.age<-substr(rodat$V1, 41, 42)
dat$pullus.age.acc<-substr(rodat$V1, 43, 43)
dat$day<-as.numeric(substr(rodat$V1, 44, 45))
dat$month<-as.numeric(substr(rodat$V1, 46, 47))
dat$year<-as.numeric(substr(rodat$V1, 48, 51))
dat$date.acc<-substr(rodat$V1, 52, 52)
dat$time<-substr(rodat$V1, 53, 56)
dat$place.code<-substr(rodat$V1, 57, 60)
dat$country<-substr(rodat$V1, 57, 58)
dat$region<-substr(rodat$V1, 59, 60)

latitude<-substr(rodat$V1, 61, 67)
latitude.grad<-as.numeric(substr(latitude, 2,3))
latitude.min<-as.numeric(substr(latitude, 4,5))
latitude.sec<-as.numeric(substr(latitude, 6,7))
latitude.sign <- substr(latitude, 1, 1)
longitude<-substr(rodat$V1, 68, 75)
longitude.grad<-as.numeric(substr(longitude, 2,4))
longitude.min<-as.numeric(substr(longitude, 5,6))
longitude.sec<-as.numeric(substr(longitude, 7,8))
longitude.sign <- substr(longitude, 1, 1)
dat$lat <- as.numeric(paste0(latitude.sign, 1)) * decimal.coord(latitude.grad+latitude.min/100)
dat$lon <- as.numeric(paste0(longitude.sign, 1)) * decimal.coord(longitude.grad+longitude.min/100)
dat$coord.acc<-substr(rodat$V1, 76, 76)


suppressWarnings(dat$condition<-as.numeric(substr(rodat$V1, 77, 77)))
suppressWarnings(dat$circumstances<-as.numeric(substr(rodat$V1, 78, 79)))
suppressWarnings(dat$circumstances.presumed<-as.numeric(substr(rodat$V1, 80, 80)))
dat$euring.codeid<-substr(rodat$V1, 81, 81)
suppressWarnings(dat$distance<-as.numeric(substr(rodat$V1, 82, 86)))
suppressWarnings(dat$direction<-as.numeric(substr(rodat$V1, 87, 89)))
suppressWarnings(dat$time.elapsed<-as.numeric(substr(rodat$V1, 90, 94)))
return(dat)
options(warn= as.numeric(warnorig))

}


