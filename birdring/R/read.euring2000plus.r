# read EURING data format EURING 2000+ into R
# Author: Fraenzi Korner-Nievergelt
#
# Reference: du Feu et al. 2012: The EURING Exchange Code 2000+. www.euring.org
#-------------------------------------------------------------------------------
read.euring2000plus <- function(filename){
#-------------------------------------------------------------------------------
# filename: name of the psv-file obtained from EURING
#-------------------------------------------------------------------------------
options(encoding="latin1")

rodat <- read.table(filename, colClasses = "character", sep="|")
dat <- data.frame(scheme=rodat$V1, id.method=rodat$V2, ring=rodat$V3)
dat$ring.verif <- rodat$V4
dat$metal.ring.info <- rodat$V5
dat$marks.info <- rodat$V6
dat$spec.byringer <- rodat$V7
dat$spec.byscheme <- rodat$V8
dat$manipulated <- rodat$V9
dat$moved <- rodat$V10
dat$catching.method <- rodat$V11
dat$catching.lures <- rodat$V12
dat$sex.byringer <- rodat$V13
dat$sex.byscheme <- rodat$V14
dat$age.byringer <- rodat$V15
dat$age.byscheme<-rodat$V16
dat$status <- rodat$V17

suppressWarnings(dat$broodsize <- as.numeric(rodat$V18))
if(sum(is.na(dat$broodsize)>0)) {
    navalues <- unique(rodat$V18[is.na(dat$broodsize)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing broodsize, the following entries were transformed into missing values (NA):\n",
navalues, "\n\n")
}

dat$pullus.age <- rodat$V19
dat$pullus.age.acc <- rodat$V20
dat$day <- as.numeric(substr(rodat$V21, 1, 2))
dat$month<-as.numeric(substr(rodat$V21, 3, 4))
dat$year<-as.numeric(substr(rodat$V21, 5, 8))
dat$date.acc <- as.numeric(rodat$V22)
dat$time<-rodat$V23
dat$place.code<- rodat$V24
dat$country<- substring(rodat$V24, 1,2)
dat$region<- substring(rodat$V24, 3,4)

latitude<-substring(rodat$V25, 1,7)
latitude.grad<-as.numeric(substr(latitude, 2,3))
latitude.min<-as.numeric(substr(latitude, 4,5))
latitude.sec<-as.numeric(substr(latitude, 6,7))
latitude.sign <- substr(latitude, 1, 1)
longitude<-substr(rodat$V25, 8, 15)
longitude.grad<-as.numeric(substr(longitude, 2,4))
longitude.min<-as.numeric(substr(longitude, 5,6))
longitude.sec<-as.numeric(substr(longitude, 7,8))
longitude.sign <- substr(longitude, 1, 1)
dat$lat <- as.numeric(paste0(latitude.sign, 1)) * decimal.coord(latitude.grad+latitude.min/100)
dat$lon <- as.numeric(paste0(longitude.sign, 1)) * decimal.coord(longitude.grad+longitude.min/100)

suppressWarnings(dat$coord.acc <- as.numeric(rodat$V26))
if(sum(is.na(dat$coord.acc)>0)) {
    navalues <- unique(rodat$V26[is.na(dat$coord.acc)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing coord.acc, the following entries were transformed into missing values (NA):\n",
navalues, "\n\n")
}

suppressWarnings(dat$condition <- as.numeric(rodat$V27))
if(sum(is.na(dat$condition)>0)){
    navalues <- unique(rodat$V27[is.na(dat$condition)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing condition, the following entries were transformed into missing values (NA):\n",
navalues, "\n\n")
}

suppressWarnings(dat$circumstances <- as.numeric(rodat$V28))
if(sum(is.na(dat$circumstances)>0)) {
    navalues <- unique(rodat$V28[is.na(dat$circumstances)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing circumstances, the following entries were transformed into missing values (NA):\n",
navalues, "\n\n")
}

suppressWarnings(dat$circumstances.presumed <- as.numeric(rodat$V29))
if(sum(is.na(dat$circumstances.presumed)>0)){
    navalues <- unique(rodat$V29[is.na(dat$circumstances.presumed)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing circumstances.presumed, the following entries were transformed into missing values (NA):\n",
navalues, "\n\n")
}

suppressWarnings(dat$euring.codeid <- as.numeric(rodat$V30))
if(sum(is.na(dat$euring.codeid)>0)){
    navalues <- unique(rodat$V30[is.na(dat$euring.codeid)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing euring.codeid, the following entries were transformed into missing values (NA):\n",
navalues, "\n\n")
}

suppressWarnings(dat$distance<-as.numeric(rodat$V31))
if(sum(is.na(dat$distance)>0)){
    navalues <- unique(rodat$V31[is.na(dat$distance)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))   cat("While importing distance, the following entries were transformed into missing values (NA):\n", 
    navalues, "\n\n")
}

suppressWarnings(dat$direction<-as.numeric(rodat$V32))
if(sum(is.na(dat$direction)>0)){
    navalues <- unique(rodat$V32[is.na(dat$direction)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing direction, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

suppressWarnings(dat$time.elapsed<-as.numeric(rodat$V33))
if(sum(is.na(dat$time.elapsed)>0)){
    navalues <- unique(rodat$V33[is.na(dat$time.elapsed)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))   cat("While importing time.elapsed, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

suppressWarnings(dat$wing.length<-as.numeric(rodat$V34))
if(sum(is.na(dat$wing.length)>0)){
    navalues <- unique(rodat$V34[is.na(dat$wing.length)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing wing.length, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

suppressWarnings(dat$third.primary <- as.numeric(rodat$V35))
if(sum(is.na(dat$third.primary)>0)){
    navalues <- unique(rodat$V35[is.na(dat$third.primary)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing third.primary, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$state.of.wing.point <- rodat$V36
suppressWarnings(dat$mass <- as.numeric(rodat$V37))
if(sum(is.na(dat$mass)>0)){
    navalues <- unique(rodat$V37[is.na(dat$mass)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing mass, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$moult <- rodat$V38
dat$moult[dat$moult==""] <- NA
dat$plumage.code <- rodat$V39
dat$plumage.code[dat$plumage.code==""] <- NA

suppressWarnings(dat$hind.claw <- as.numeric(rodat$V40))
if(sum(is.na(dat$hind.claw)>0)) {
    navalues <- unique(rodat$V40[is.na(dat$hind.claw)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))cat("While importing hind.claw, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

suppressWarnings(dat$bill.length <- as.numeric(rodat$V41))
if(sum(is.na(dat$bill.length)>0)) {
    navalues <- unique(rodat$V41[is.na(dat$bill.length)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing bill.length, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$bill.method <- rodat$V42
dat$bill.method[dat$bill.method==""] <- NA

suppressWarnings(dat$total.head.length <- as.numeric(rodat$V43))
if(sum(is.na(dat$total.head.length)>0)){
    navalues <- unique(rodat$V43[is.na(dat$total.head.length)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing total.head.length, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

suppressWarnings(dat$tarsus <- as.numeric(rodat$V44))
if(sum(is.na(dat$tarsus)>0)){
    navalues <- unique(rodat$V44[is.na(dat$tarsus)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing tarsus, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$tarsus.method <- rodat$V45
dat$tarsus.method[dat$tarsus.method==""] <- NA
suppressWarnings(dat$tail.length <- as.numeric(rodat$V46))
if(sum(is.na(dat$tail.length)>0)){
    navalues <- unique(rodat$V46[is.na(dat$tail.length)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing tail.length, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

suppressWarnings(dat$tail.difference <- as.numeric(rodat$V47))
if(sum(is.na(dat$tail.difference)>0)) {
    navalues <- unique(rodat$V47[is.na(dat$tail.difference)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing tail.difference, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

suppressWarnings(dat$fat.score <- as.numeric(rodat$V48))
if(sum(is.na(dat$fat.score)>0))  {
    navalues <- unique(rodat$V48[is.na(dat$fat.score)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing fat.score, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}


dat$fat.score.method <- rodat$V49
dat$fat.score.method[dat$fat.score.method==""] <- NA
suppressWarnings(dat$pectoral.muscle <- as.numeric(rodat$V50))
if(sum(is.na(dat$pectoral.muscle)>0)) {
    navalues <- unique(rodat$V50[is.na(dat$pectoral.muscle)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0)) cat("While importing pectoral.muscle, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$brood.patch <- rodat$V51
suppressWarnings(dat$primary.score <- as.numeric(rodat$V52))
if(sum(is.na(dat$primary.score)>0)){
    navalues <- unique(rodat$V52[is.na(dat$primary.score)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing primary.score, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$primary.moult <- rodat$V53
dat$primary.moult[dat$primary.moult==""] <- NA
suppressWarnings(dat$old.greater.coverts <- as.numeric(rodat$V54))
if(sum(is.na(dat$old.greater.coverts )>0)){
    navalues <- unique(rodat$V54[is.na(dat$old.greater.coverts )]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing old.greater.coverts , the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$alula <- rodat$V55
dat$alula[dat$alula==""] <- NA
suppressWarnings(dat$carpal.covert <- as.numeric(rodat$V56))
if(sum(is.na(dat$carpal.covert)>0)){
    navalues <- unique(rodat$V56[is.na(dat$carpal.covert)]); navalues <- navalues[navalues!=""]
    if(length(navalues>0))  cat("While importing carpal.covert, the following entries were transformed into missing values (NA):\n",
    navalues, "\n\n")
}

dat$sexing.method <- rodat$V57
dat$sexing.method[dat$sexing.method==""] <- NA
if(ncol(rodat)>57){
dat$place.name <- rodat$V58
dat$remarks <- rodat$V59
dat$reference <- rodat$V60
}

return(dat)

}

