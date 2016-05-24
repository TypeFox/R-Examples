rm(list=ls())
namstates <- c("H","A","B","C")
nsample <- 22
id <- 1:nsample
# Dates of birth: randomly distributed in 1991
#   = random number from uniform distribution      
birth <- as.Date("1991-01-01")+runif(22,min=0,max=365)  
# Observation window
entry = c("02/01/2007", "17/01/2007", "18/01/2007","22/01/2007","10/02/2007",
"30/01/2007","04/04/2007","29/04/2007","18/05/2007","20/05/2007","15/05/2007",
"05/02/2007","05/02/2007","06/02/2007","26/02/2007","10/03/2007",
"11/03/2007","28/03/2007","15/03/2007","13/04/2007","04/04/2007","25/04/2007")
interview = c("25/05/2007", "17/05/2007", "10/05/2007","13/05/2007", "23/05/2007",
"15/05/2007","06/05/2007", "27/05/2007", "29/05/2007","31/05/2007", "18/05/2007",
"19/05/2007","10/05/2007", "28/05/2007", "22/05/2007","25/05/2007", "12/05/2007",
"29/05/2007","10/05/2007", "20/05/2007", "11/05/2007","31/05/2007")
# Covariates
cov=rep("X",length(id))
# Create a data frame of dates of transition
trans1 = c("11/02/2007", "04/05/2007", "NA","28/02/2007", "17/05/2007",
"12/02/2007","NA", "NA", "NA","NA", "NA",
"25/02/2007","18/04/2007", "18/05/2007", "NA","NA", "08/05/2007",
"NA","23/03/2007", "NA", "09/05/2007","16/05/2007")
trans2 =c ("23/03/2007", "NA", "NA","10/04/2007", "NA",
"05/03/2007","NA", "NA", "NA","NA", "NA",
"01/04/2007","30/04/2007", "20/05/2007", "NA","NA", "NA",
"NA","08/04/2007", "NA", "NA","20/05/2007")
trans3 =c("05/05/2007", "NA", "NA","10/05/2007", "NA",
"17/04/2007","NA", "NA", "NA","NA", "NA",
"02/05/2007","NA", "NA", "NA","NA", "NA",
"NA","20/04/2007", "NA", "NA","26/05/2007")
h=data.frame (A=trans1,B=trans2,C=trans3,stringsAsFactors=FALSE)
tp <- data.frame (ID=id,Born=as.character(format(birth,"%d/%m/%Y")),Start=entry,Stop=interview,A=trans1,B=trans2,C=trans3,stringsAsFactors=FALSE)
# Convert character object d to object of class 'Date' (days elapsed since 1-1-1970 (Julian dates))
d<- apply(h,2,function(x) y=as.Date(x,"%d/%m/%Y"))
d <- data.frame(d)  #  d is numeric
dimnames(d) <- dimnames(h)
# nsample = sample size
nsample <- nrow(d)
# Sort and rearrange the dates and obtain state sequence (path)
f <- Sequences.ind.0(d=d,namstates=namstates,absorb="C")
# Convert Julian dates in calendar (Gregorian) dates
dates <- data.frame (f$d)
for (i in 1:3)
 {dates[,i] <- as.Date(dates[,i],origin="1970-01-01") 
 }
colnames(dates) <- paste ("Tr",1:(length(namstates)-1),sep="")
# State sequence
path <- as.character(f$path)
maxns <- max (nchar(path)) 

RS <- data.frame (ID=id,born=birth,start=as.Date(entry,"%d/%m/%Y"),end=as.Date(interview,"%d/%m/%Y"),cov=cov,path=as.character(path),dates[,1:(maxns-1)],stringsAsFactors=FALSE)
attr(RS,"format.date") <- "%Y-%m-%d"
attr(RS,"format.born") <- "%Y-%m-%d"
attr(RS,"param") <- Parameters (RS)
#class (RS) <- c("Biograph","data.frame")
zz8 <- "/Users/frans/Documents/R/0 0 MAC/Package/Biograph.TEST/Chapters/AnnexA/simple 2/"
setwd(zz8)
save (RS,file=paste(zz8,"simple2.RData",sep=""))
# ==============  Biograph object created  ====================

